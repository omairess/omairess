# ==============================================================================
# tests/testthat/test-time-origin.R
#
# User report: with Time origin = "first observation" the fit plot started at
# 00:00 for a recording that began at 08:00, and with "midnight" it looked
# right -- so the toggle appeared to work backwards.
#
# It was not backwards. The MODEL was correct either way; the DISPLAY was not.
# time_origin = "first_observation" re-anchors the fit so mod$time_vec holds
# 0, 1, ... 31 for a recording that ran 08:00 -> 15:00 two days later. Anything
# that labels that axis as a clock reads t = 0 as midnight, so the axis sat
# eight hours early. Under "midnight" the shift is zero and model time IS clock
# time, which is why only one of the two settings looked wrong.
#
# The fix is one addition at every DISPLAY boundary -- tick text, hover, and the
# angle a polar point is drawn at -- and never inside a fit. These tests pin
# both halves: the model keeps its own origin, and everything shown to a reader
# is converted back.
# ==============================================================================

`%||%` <- function(a, b) if (is.null(a)) b else a
app_dir <- if (dir.exists("server")) "." else if (dir.exists("../../server")) "../.." else "FCK"
source(file.path(app_dir, "server/08_helpers_cosinor.R"))

# the user's actual axis: 19 points, 08:00 -> 15:00 two days later
clock_times <- c(8, 9, 10, 11, 12, 14, 16, 18, 20, 21, 22, 23, 0, 2, 4, 6, 8, 10, 15)
linear_times <- c(8, 9, 10, 11, 12, 14, 16, 18, 20, 21, 22, 23, 24, 26, 28, 30, 32, 34, 39)

mk_mod <- function(origin) {
  shift <- if (identical(origin, "first_observation")) min(linear_times) else 0
  list(period = 24, time_origin = origin, origin_shift = shift,
       time_vec = linear_times - shift,      # what the fit runs on
       time_vec_clock = linear_times)        # what a reader should see
}

test_that("the two origins differ by exactly the start hour", {
  mid <- mk_mod("midnight"); fst <- mk_mod("first_observation")
  expect_equal(mid$origin_shift, 0)
  expect_equal(fst$origin_shift, 8)
  expect_equal(mid$time_vec, linear_times)
  expect_equal(fst$time_vec, linear_times - 8)
  expect_equal(fst$time_vec[1], 0)           # the model starts at zero
  expect_equal(mid$time_vec[1], 8)
})

test_that("model time converts back to the SAME clock time under either origin", {
  # this is the invariant the bug broke
  for (o in c("midnight", "first_observation")) {
    mod <- mk_mod(o)
    expect_equal(fck_model_to_clock(mod$time_vec, mod), linear_times, info = o)
  }
})

test_that("the first observation is 08:00 on the clock whichever origin is set", {
  for (o in c("midnight", "first_observation")) {
    mod <- mk_mod(o)
    lab <- fck_clock_label(fck_model_to_clock(mod$time_vec, mod), 24)
    expect_equal(lab[1], "08:00", info = o)
  }
})

test_that("the unconverted axis is what produced the reported symptom", {
  fst <- mk_mod("first_observation")
  # the bug: labelling model time directly
  wrong <- fck_clock_label(fst$time_vec, 24)
  expect_equal(wrong[1], "00:00")            # exactly what the user saw
  # the fix
  right <- fck_clock_label(fck_model_to_clock(fst$time_vec, fst), 24)
  expect_equal(right[1], "08:00")
  expect_false(identical(wrong[1], right[1]))
  # and under "midnight" the two coincide, which is why only one looked wrong
  mid <- mk_mod("midnight")
  expect_equal(fck_clock_label(mid$time_vec, 24),
               fck_clock_label(fck_model_to_clock(mid$time_vec, mid), 24))
})

test_that("ticks land on the model axis but read in clock time", {
  fst <- mk_mod("first_observation")
  shift <- fst$origin_shift
  tk <- fck_clock_ticks(range(fst$time_vec) + shift, 24)
  vals_on_model_axis <- tk$vals - shift
  # positions stay inside the model's own range, so they plot where the data are
  expect_true(all(vals_on_model_axis >= min(fst$time_vec) - 1e-9 &
                  vals_on_model_axis <= max(fst$time_vec) + 1e-9))
  # and the labels are clock times covering the recording. The first tick sits
  # on the first multiple of the step at or after the start, so for an 08:00
  # start with a 3 h step it is 09:00 -- the tick grid is regular, it does not
  # begin at the data.
  expect_true(all(nchar(tk$text) == 5))
  expect_gte(tk$vals[1], min(linear_times))
  expect_lt(tk$vals[1] - min(linear_times), tk$step)
  expect_equal(tk$text[1], fck_clock_label(tk$vals[1], 24))
})

test_that("the clock range covers the real recording under either origin", {
  for (o in c("midnight", "first_observation")) {
    mod <- mk_mod(o)
    expect_equal(fck_model_clock_range(mod), range(linear_times), info = o)
  }
})

test_that("midnight boundaries are clock events, not model events", {
  fst <- mk_mod("first_observation")
  shift <- fst$origin_shift
  cr <- range(fst$time_vec) + shift              # 8 .. 39
  bnds_clock <- seq(ceiling(cr[1] / 24) * 24, cr[2], by = 24)
  bnds_clock <- bnds_clock[bnds_clock > cr[1] & bnds_clock < cr[2]]
  expect_equal(bnds_clock, 24)                   # one midnight, at clock 24:00
  # drawn on the model axis it sits at 16, not 24
  expect_equal(bnds_clock - shift, 16)
  # naively using model time would put a "midnight" at model 24, which is
  # clock 32:00 -- 08:00 the next morning, not midnight at all
  naive <- seq(24, max(fst$time_vec), by = 24)
  expect_equal(naive, 24)
  expect_equal(fck_clock_label(naive + shift, 24), "08:00")
  expect_false(identical(naive, bnds_clock - shift))
})

test_that("a polar dial places a point at its CLOCK angle", {
  fst <- mk_mod("first_observation")
  shift <- fst$origin_shift
  t_abs <- fst$time_vec                        # model time
  right <- (t_abs + shift) %% 24
  wrong <- t_abs %% 24
  # the first observation belongs at 08:00 on the dial
  expect_equal(right[1], 8)
  expect_equal(wrong[1], 0)                    # the bug: eight hours early
  expect_equal(right, linear_times %% 24)
})

test_that("laps are counted in clock days, not model days", {
  fst <- mk_mod("first_observation")
  shift <- fst$origin_shift
  t_abs <- seq(min(fst$time_vec), max(fst$time_vec), length.out = 400)
  lap_right <- floor((t_abs + shift) / 24)
  lap_wrong <- floor(t_abs / 24)
  # both give two laps here, so the COUNT does not reveal the error -- the
  # BOUNDARY does. On clock time the split falls at midnight; on model time it
  # falls at 08:00, cutting the recording in the wrong place.
  expect_equal(length(unique(lap_right)), 2)
  split_right <- t_abs[which(diff(lap_right) != 0) + 1]
  split_wrong <- t_abs[which(diff(lap_wrong) != 0) + 1]
  expect_equal(fck_clock_label(split_right + shift, 24), "00:00", tolerance = 0)
  expect_equal(fck_clock_label(split_wrong + shift, 24), "08:00")
  expect_false(isTRUE(all.equal(split_right, split_wrong)))
})

test_that("nothing in the conversion touches the fitted values", {
  # the trend is what forces linear time, and it must be evaluated on MODEL
  # time; converting for display must not change it
  fst <- mk_mod("first_observation")
  S_model <- 1 - exp(-(fst$time_vec - min(fst$time_vec)) / 18)
  S_clock <- 1 - exp(-(fck_model_to_clock(fst$time_vec, fst) -
                       min(fck_model_to_clock(fst$time_vec, fst))) / 18)
  expect_equal(S_model, S_clock)              # elapsed time is origin-free
  expect_true(all(diff(S_model) > 0))         # and monotone, as it must be
})

test_that("a model with no origin_shift field behaves as unshifted", {
  # backwards compatibility with a saved session from before the option existed
  mod <- list(period = 24, time_vec = linear_times)
  expect_equal(fck_model_to_clock(mod$time_vec, mod), linear_times)
  expect_equal(fck_model_clock_range(mod), range(linear_times))
})
