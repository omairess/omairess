# ==============================================================================
# tests/testthat/test-traj-acrophase-display.R
#
# Reported from a real run: the H1 acrophase column in tab 6 was EMPTY, and the
# acrophases that were visible elsewhere looked wrong against the plotted curves.
#
# Two separate things, one a bug and one not.
#
# THE BUG. The renderer read dance_acrophase_clock(...)$first. That function
# returns hours / all_hours / elapsed / effective_period / harmonic /
# clock_origin -- there is no $first; the name came from the COMMENT on the first
# element, "first maximum on the clock". R answers a missing name with NULL, and
# dance_clock_label(NULL) is character(0), and tags$td(character(0)) is an empty
# cell. No error, no NA, no warning: a silently blank column. A test that only
# asked whether the table rendered could not see it, so the check below asks the
# parser which field every call site actually reads.
#
# NOT A BUG. An H1 acrophase is the maximum of the first harmonic alone. The
# plotted trajectory is level + trend + H1 + H2, and H2 has its own maximum, so
# the drawn curve peaks somewhere else -- measured below at 1.3 to 1.9 h earlier
# for a second harmonic around 40% of the first. Both numbers are right; they
# are answers to different questions, and the app now shows both.
# ==============================================================================
.libPaths(c("~/Rlib", .libPaths()))
app_dir <- if (dir.exists("server")) "." else if (dir.exists("../../server")) "../.." else "DANCE"
`%||%` <- function(a, b) if (is.null(a)) b else a

e <- new.env()
for (f in c("08_helpers_cosinor", "08c_helpers_circstat", "08d_helpers_traj",
            "08e_helpers_trajfit", "08f_helpers_trajinf", "08g_helpers_trajband"))
  sys.source(file.path(app_dir, "server", paste0(f, ".R")), envir = e)

test_that("the failure mode is a BLANK cell, which is why it went unseen", {
  a <- e$dance_acrophase_clock(hours = 3, clock_origin = 8)
  expect_null(a$first)                       # the name the renderer used
  expect_equal(a$hours, 11)                  # the name that exists
  # NULL travels all the way to an empty cell without ever erroring
  expect_length(e$dance_clock_label(a$first, 24), 0L)
  expect_equal(e$dance_clock_label(a$hours, 24), "11:00")
})

test_that("every call site reads a field dance_acrophase_clock actually returns", {
  # Asked of the PARSER, not of the text: a regex cannot tell the $period in
  # `ff$spec$period` (an argument) from the $hours on the result, and pinning the
  # spelling of one call site is what let the original mistake through.
  valid <- names(e$dance_acrophase_clock(hours = 1))
  expect_true(all(c("hours", "all_hours", "elapsed") %in% valid))
  found <- list()
  walk <- function(x) {
    if (is.call(x)) {
      if (length(x) == 3L && identical(x[[1]], as.name("$")) &&
          is.call(x[[2]]) && identical(x[[2]][[1]], as.name("dance_acrophase_clock")))
        found[[length(found) + 1L]] <<- as.character(x[[3]])
      for (i in seq_along(x)) if (!is.null(x[[i]])) try(walk(x[[i]]), silent = TRUE)
    }
  }
  files <- list.files(file.path(app_dir, "server"), pattern = "[.]R$", full.names = TRUE)
  for (f in files) for (ex in parse(f)) walk(ex)
  expect_gt(length(found), 0L)               # the accessor exists somewhere
  expect_true(all(unlist(found) %in% valid),
              info = paste("invalid field(s):",
                           paste(setdiff(unlist(found), valid), collapse = ", ")))
})

# ---------------------------------------------------------------- behaviour --
mk <- function(nper = 12, seed = 20) {
  set.seed(seed)
  P <- 24; tp <- seq(8, 32, length.out = 13)
  truth <- c(A = 2.0, B = 1.0)               # H1 peak, clock hours
  grp <- rep(names(truth), each = nper)
  Y <- t(sapply(seq_along(grp), function(i) {
    g <- grp[i]
    50 + rnorm(1, 0, 6) +
      16 * cos(2 * pi * (tp - truth[[g]]) / P) +
      7  * cos(2 * pi * 2 * (tp - (truth[[g]] - 3)) / P) +
      0.20 * (tp - 8) + rnorm(length(tp), 0, 4)
  }))
  sp <- e$dance_traj_spec(e$dance_traj_long(
    Y, tp, sprintf("S%03d", seq_along(grp)), list(Group = grp)), 24, 2, "linear")
  ff <- e$dance_traj_fit(sp)
  ff$clock_origin <- sp$t0                   # raw clock times in, so t0 is the origin
  list(fit = ff, truth = truth, period = P)
}

test_that("the H1 acrophase recovers a PLANTED clock time", {
  skip_if_not_installed("lme4"); skip_if_not_installed("emmeans")
  m <- mk()
  skip_if_not(isTRUE(m$fit$ok))
  co <- e$dance_traj_cell_coefs(m$fit, 1)
  skip_if_not(isTRUE(co$ok))
  ap <- e$dance_traj_amp_phase_ci(co, method = "joint", n_draw = 6000)
  skip_if_not(isTRUE(ap$ok))
  for (i in seq_len(nrow(ap$table))) {
    r <- ap$table[i, ]
    got <- e$dance_acrophase_clock(hours = r$acrophase_time, period = m$period,
                                   harmonic = 1, clock_origin = m$fit$clock_origin)$hours
    want <- m$truth[[r$cell]]
    gap <- abs(((got - want + 12) %% 24) - 12)
    expect_lt(gap, max(0.5, r$acrophase_arc_time))   # inside its own arc
    expect_match(e$dance_clock_label(got, m$period, show_day = FALSE),
                 "^[0-2][0-9]:[0-5][0-9]$")          # and rendered as hh:mm
  }
})

test_that("the drawn curve peaks EARLIER than H1, because H2 is in the model", {
  skip_if_not_installed("lme4"); skip_if_not_installed("emmeans")
  m <- mk()
  skip_if_not(isTRUE(m$fit$ok))
  pk <- e$dance_traj_curve_peaks(m$fit, "full")
  expect_true(isTRUE(pk$ok))
  expect_equal(nrow(pk$table), 2L)
  expect_true(all(pk$table$peak_interior))   # a real turning point, not the edge

  # independent oracle: evaluate the fitted curve and take the argmax directly
  bv <- e$dance_traj_beta(m$fit); sp <- m$fit$spec
  grid <- e$dance_traj_cell_grid(sp)
  tt <- seq(sp$t0, sp$t0 + sp$period, length.out = 2001)
  for (i in seq_len(nrow(grid))) {
    X <- e$dance_traj_design_rows(m$fit, grid[i, , drop = FALSE], tt)
    keep <- intersect(colnames(X), names(bv$beta))
    y <- as.numeric(X[, keep, drop = FALSE] %*% bv$beta[keep])
    want <- tt[which.max(y)]
    got <- pk$table$peak_t[match(grid$.cell[i], pk$table$cell)]
    expect_lt(abs(got - want), 0.05)         # within the grid resolution
  }

  co <- e$dance_traj_cell_coefs(m$fit, 1)
  ap <- e$dance_traj_amp_phase_ci(co, method = "joint", n_draw = 4000)
  gaps <- vapply(seq_len(nrow(ap$table)), function(i) {
    h1 <- e$dance_acrophase_clock(hours = ap$table$acrophase_time[i], period = m$period,
                                  harmonic = 1, clock_origin = m$fit$clock_origin)$hours
    pkc <- pk$table$peak_t[match(ap$table$cell[i], pk$table$cell)] %% m$period
    ((h1 - pkc + 12) %% 24) - 12
  }, numeric(1))
  # H1 LEADS the drawn peak here, by more than any acrophase interval
  expect_true(all(gaps > 0.5),
              info = sprintf("gaps: %s", paste(round(gaps, 2), collapse = ", ")))
  expect_true(all(gaps < 4))
})

test_that("the peak search is confined to one period and says its resolution", {
  skip_if_not_installed("lme4")
  m <- mk()
  skip_if_not(isTRUE(m$fit$ok))
  pk <- e$dance_traj_curve_peaks(m$fit, "full")
  expect_lte(diff(pk$window), m$fit$spec$period + 1e-9)
  expect_true(all(pk$table$peak_t >= pk$window[1] & pk$table$peak_t <= pk$window[2]))
  expect_lt(pk$resolution_min, 2)            # sub-two-minute grid
})

# ==============================================================================
# The peak row shipped with the clock offset LEFT OUT: the table converted
# model-elapsed hours with `peak_t %% period`, so on a study starting at 08:00 a
# peak at 08:00 (+1d) was reported as 00:00. The previous test checked peak_t --
# the elapsed value, which was right -- and never the number a reader sees.
# The conversion now lives in the helper, and these check it there.
# ==============================================================================

# Two fixtures, because the formula has two parts and ONE fixture can only pin
# one of them. The app builds its model on an axis that already starts at zero
# and carries the clock offset separately (t0 = 0, origin = 8), so a conversion
# that forgets t0 survives every run inside the app. A fit on raw clock times
# has t0 = 8 AND origin = 8, where the two cancel, so a conversion that forgets
# the origin survives that one. Only both together fix (t - t0) + origin.
mk_elapsed <- function(nper = 12, seed = 20) {
  m <- mk(nper, seed)
  set.seed(seed)
  P <- 24; tp <- seq(0, 24, length.out = 13)      # ELAPSED, as the app builds it
  truth <- c(A = 2.0, B = 1.0)
  grp <- rep(names(truth), each = nper)
  Y <- t(sapply(seq_along(grp), function(i) {
    g <- grp[i]
    50 + rnorm(1, 0, 6) + 16 * cos(2*pi*((tp + 8) - truth[[g]])/P) +
      7 * cos(2*pi*2*((tp + 8) - (truth[[g]] - 3))/P) + 0.20 * tp + rnorm(length(tp), 0, 4)
  }))
  sp <- e$dance_traj_spec(e$dance_traj_long(
    Y, tp, sprintf("S%03d", seq_along(grp)), list(Group = grp)), 24, 2, "linear")
  ff <- e$dance_traj_fit(sp)
  ff$clock_origin <- 8                            # the offset back to the clock
  list(fit = ff, truth = truth, period = P)
}

test_that("the peak carries its own clock time, measured from the BASIS origin", {
  skip_if_not_installed("lme4")
  for (nm in c("raw clock axis", "elapsed axis, as the app builds it")) {
    m <- if (identical(nm, "raw clock axis")) mk() else mk_elapsed()
    skip_if_not(isTRUE(m$fit$ok))
    sp <- m$fit$spec; P <- sp$period; co <- m$fit$clock_origin
    pk <- e$dance_traj_curve_peaks(m$fit, "full")
    expect_true(isTRUE(pk$ok), info = nm)
    for (i in seq_len(nrow(pk$table))) {
      r <- pk$table[i, ]
      expect_equal(r$peak_clock, (r$peak_t - sp$t0 + co) %% P, info = nm)
      expect_equal(r$trough_clock, (r$trough_t - sp$t0 + co) %% P, info = nm)
    }
    expect_equal(pk$window_clock, (pk$window - sp$t0 + co) %% P, info = nm)
    if (!identical(nm, "raw clock axis")) {
      # on the app's axis the origin is what the old code dropped, and dropping
      # it must give a DIFFERENT answer or this proves nothing
      expect_equal(sp$t0, 0)
      expect_false(isTRUE(all.equal(pk$table$peak_clock, pk$table$peak_t %% P)))
    } else {
      # here t0 and the origin cancel, which is why the raw-clock fixture alone
      # could never have caught the missing origin
      expect_equal(sp$t0, co)
    }
  }
})

test_that("with one harmonic and no trend, the curve peak IS the acrophase", {
  skip_if_not_installed("lme4"); skip_if_not_installed("emmeans")
  # the one configuration where the two are the same number, so the peak search
  # can be checked against an independent quantity rather than against itself
  set.seed(44)
  tp <- seq(8, 32, length.out = 25)
  grp <- rep(c("A", "B"), each = 14)
  Y <- t(sapply(seq_along(grp), function(i)
    40 + rnorm(1, 0, 4) + 12 * cos(2 * pi * (tp - 2) / 24) + rnorm(length(tp), 0, 2)))
  sp <- e$dance_traj_spec(e$dance_traj_long(
    Y, tp, sprintf("S%03d", seq_along(grp)), list(Group = grp)), 24, 1, "none")
  ff <- e$dance_traj_fit(sp); skip_if_not(isTRUE(ff$ok))
  ff$clock_origin <- sp$t0
  pk <- e$dance_traj_curve_peaks(ff, "full")
  co <- e$dance_traj_cell_coefs(ff, 1); skip_if_not(isTRUE(co$ok))
  ap <- e$dance_traj_amp_phase_ci(co, method = "joint", n_draw = 4000)
  for (i in seq_len(nrow(ap$table))) {
    acro <- e$dance_acrophase_clock(hours = ap$table$acrophase_time[i], period = 24,
                                    harmonic = 1, clock_origin = ff$clock_origin)$hours
    peak <- pk$table$peak_clock[match(ap$table$cell[i], pk$table$cell)]
    gap <- abs(((acro - peak + 12) %% 24) - 12)
    expect_lt(gap, 0.1)                      # within the search grid
    expect_lt(abs(((acro - 2 + 12) %% 24) - 12), 0.5)   # and both near the truth
  }
})

test_that("nothing renders a peak by converting the elapsed value itself", {
  # the bug was a call site, not the helper, so the call sites are what is checked
  src <- paste(vapply(parse(file.path(app_dir, "server/72_harmonic.R")),
                      function(x) paste(deparse(x), collapse = " "), character(1)),
               collapse = " ")
  src <- gsub("[[:space:]]+", " ", src)
  expect_false(grepl("peak_t %% ", src, fixed = TRUE))
  expect_true(grepl("peak_clock", src, fixed = TRUE))
})

test_that("the acrophase interval endpoints span exactly the reported arc", {
  skip_if_not_installed("lme4"); skip_if_not_installed("emmeans")
  m <- mk(); skip_if_not(isTRUE(m$fit$ok))
  co <- e$dance_traj_cell_coefs(m$fit, 1); skip_if_not(isTRUE(co$ok))
  ap <- e$dance_traj_amp_phase_ci(co, method = "joint", n_draw = 6000)
  P <- m$period
  for (i in seq_len(nrow(ap$table))) {
    r <- ap$table[i, ]
    if (!isTRUE(r$phase_defined)) next
    # the arc is the width of the interval the endpoints describe -- the table
    # shows both, and they have to be the same statement
    expect_equal((r$acrophase_hi - r$acrophase_lo) %% P, r$acrophase_arc_time,
                 tolerance = 1e-8)
    expect_gt(r$acrophase_arc_time, 0)
    expect_lt(r$acrophase_arc_time, P)
    # the estimate lies inside its own interval, the short way round
    d1 <- ((r$acrophase_time - r$acrophase_lo) %% P)
    expect_lte(d1, r$acrophase_arc_time + 1e-8)
  }
})

test_that("an undefined phase reports no interval rather than a misleading one", {
  skip_if_not_installed("lme4"); skip_if_not_installed("emmeans")
  # a fixture built so H2 IS weak: a second harmonic well inside the noise, so
  # its (cos, sin) region covers the origin and its angle means nothing
  set.seed(77)
  tp <- seq(0, 24, length.out = 13); nper <- 10
  grp <- rep(c("A", "B"), each = nper)
  Y <- t(sapply(seq_along(grp), function(i)
    40 + rnorm(1, 0, 5) + 14 * cos(2*pi*tp/24) + 0.15 * cos(2*pi*2*tp/24) +
      rnorm(length(tp), 0, 6)))
  sp <- e$dance_traj_spec(e$dance_traj_long(
    Y, tp, sprintf("S%03d", seq_along(grp)), list(Group = grp)), 24, 2, "none")
  ff <- e$dance_traj_fit(sp); skip_if_not(isTRUE(ff$ok))
  ff$clock_origin <- 8
  co <- e$dance_traj_cell_coefs(ff, 2)
  skip_if_not(isTRUE(co$ok))
  ap <- e$dance_traj_amp_phase_ci(co, method = "joint", n_draw = 4000)
  und <- !ap$table$phase_defined
  expect_true(any(und))          # the fixture must actually produce one
  expect_true(all(is.na(ap$table$acrophase_lo[und])))
  expect_true(all(is.na(ap$table$acrophase_hi[und])))
  expect_true(all(is.na(ap$table$acrophase_arc_time[und])))
  # the POINT estimate still exists -- it is the interval that does not
  expect_true(all(is.finite(ap$table$acrophase_time[und])))
})
