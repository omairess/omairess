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
