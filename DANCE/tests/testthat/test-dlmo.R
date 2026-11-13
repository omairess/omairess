# ==============================================================================
# tests/testthat/test-dlmo.R
#
# DLMO on profiles whose onset is known by construction: the interpolated
# crossing is exact on a linear rise, a single sample above the threshold is
# not an onset, censored profiles get no time, gaps are not filled in, two
# nights in one row are two profiles on the same clock, and the phase angle is
# bedtime minus onset. The comparisons recover a known night difference.
# ==============================================================================
.libPaths(c("~/Rlib", .libPaths()))
app_dir <- if (dir.exists("server")) "." else if (dir.exists("../../server")) "../.." else "DANCE"
`%||%` <- function(a, b) if (is.null(a)) b else a
e <- new.env()
for (f in c("08_helpers_cosinor", "08l_helpers_dlmo"))
  sys.source(file.path(app_dir, "server", paste0(f, ".R")), envir = e)

tt <- 18:26
rise <- function(on = 20, slope = 4, base = 1) base + slope * pmax(0, tt - on)

test_that("a linear rise is crossed exactly where it reaches the threshold", {
  r <- e$dance_dlmo_cross(tt, rise(), 3)
  expect_identical(r$status, "determined")
  expect_equal(r$time, 20.5)
  expect_equal(c(r$t_before, r$t_after), c(20, 21))
  # the threshold is the user's: 5 pg/mL is reached at 21
  expect_equal(e$dance_dlmo_cross(tt, rise(), 5)$time, 21)
})

test_that("one sample above the threshold is not an onset", {
  y <- rise(); y[2] <- 6
  r <- e$dance_dlmo_cross(tt, y, 3, confirm = 1)
  expect_equal(r$time, 20.5); expect_equal(r$spikes, 1L)
  # with no confirmation asked, the spike is the onset
  expect_equal(e$dance_dlmo_cross(tt, y, 3, confirm = 0)$time, 18 + 2 / 5)
})

test_that("censored profiles get no time, and say which side", {
  expect_identical(e$dance_dlmo_cross(tt, rep(10, 9), 3)$status, "before_window")
  expect_true(is.na(e$dance_dlmo_cross(tt, rep(10, 9), 3)$time))
  expect_identical(e$dance_dlmo_cross(tt, rep(1, 9), 3)$status, "no_onset")
  # a crossing at the last sample cannot be confirmed by a later one
  r <- e$dance_dlmo_cross(tt, c(rep(1, 8), 9), 3)
  expect_identical(r$status, "unconfirmed"); expect_equal(r$time, 25.25)
  expect_identical(e$dance_dlmo_cross(1, 5, 3)$status, "few_samples")
})

test_that("a missing sample is not filled in: interpolation spans the observed neighbours", {
  y <- rise(); y[4] <- NA          # 21:00 missing
  r <- e$dance_dlmo_cross(tt, y, 3)
  # between 20:00 (1) and 22:00 (9): 3 is reached a quarter of the way
  expect_equal(r$time, 20.5)
  expect_identical(r$status, "determined")
  expect_identical(e$dance_dlmo_cross(tt, y, 3, max_gap = 1)$status, "gap")
})

test_that("variable threshold: baseline mean + 2 SD of the first three samples", {
  y <- c(1, 1.2, 0.8, 1, 2, 4, 8, 12, 14)
  expect_equal(e$dance_dlmo_variable_threshold(y, 3, 2), 1.4)
  # missing samples are skipped, not counted as baseline
  expect_equal(e$dance_dlmo_variable_threshold(c(NA, y), 3, 2), 1.4)
  expect_true(is.na(e$dance_dlmo_variable_threshold(c(1, NA), 3, 2)))
  p <- list(id = "a", night = "1", group = NA, row = 1, t = tt, y = y)
  r <- e$dance_dlmo_estimate(p, "variable", list(n_base = 3, k_sd = 2))
  expect_equal(r$threshold, 1.4)
  expect_equal(r$dlmo_h, 21 + 0.4)        # 1 at 21:00, 2 at 22:00
})

test_that("values below the LLOQ are replaced by the chosen rule, and counted", {
  y <- c(0.2, 0.5, 2, 6)
  expect_equal(e$dance_dlmo_lloq(y, 1, "half")$y, c(0.5, 0.5, 2, 6))
  expect_equal(e$dance_dlmo_lloq(y, 1, "lloq")$y, c(1, 1, 2, 6))
  expect_equal(e$dance_dlmo_lloq(y, 1, "zero")$n_below, 2L)
  expect_equal(e$dance_dlmo_lloq(y)$y, y)
})

test_that("two nights in one row are two profiles read on the same clock", {
  times <- c(18:26, 42:50)
  Y <- rbind(c(rise(20), rise(21)), c(rise(19.5), rise(22)))
  pr <- e$dance_dlmo_profiles(Y, times, c("P1", "P2"))
  expect_length(pr, 4)
  expect_setequal(unique(vapply(pr, `[[`, "", "night")), c("Night 1", "Night 2"))
  tab <- e$dance_dlmo_run(pr, "fixed", list(threshold = 3))
  expect_equal(attr(tab, "anchor"), 18)
  expect_equal(tab$dlmo_h, c(20.5, 21.5, 20, 22.5))
  expect_identical(tab$dlmo_clock, c("20:30", "21:30", "20:00", "22:30"))
  # one row per night, labelled by a column
  pr2 <- e$dance_dlmo_profiles(rbind(rise(20), rise(21)), tt, c("P1", "P1"), night = c("baseline", "SD"))
  expect_identical(vapply(pr2, `[[`, "", "night"), c("baseline", "SD"))
})

test_that("onsets after midnight stay after the evening", {
  tt2 <- c(22, 23, 24, 25, 26)
  pr <- e$dance_dlmo_profiles(rbind(c(1, 1, 1, 5, 9)), tt2, "P")
  tab <- e$dance_dlmo_run(pr, "fixed", list(threshold = 3))
  expect_equal(tab$dlmo_h, 24.5)
  expect_identical(tab$dlmo_clock, "00:30")
})

test_that("phase angle = bedtime - DLMO, with the bedtime within 12 h of the onset", {
  expect_equal(e$dance_dlmo_phase_angle(21, 23), 2)
  expect_equal(e$dance_dlmo_phase_angle(25, 23.5), -1.5)
  expect_equal(e$dance_dlmo_phase_angle(21.5, 0.5), 3)
  expect_equal(e$dance_dlmo_clock_hours(c("23:30", "00:15", "23.5", "0.9791667", NA, "25:00")),
               c(23.5, 0.25, 23.5, 23.5, NA, NA), tolerance = 1e-6)
  # without a bedtime there is no phase angle
  pr <- e$dance_dlmo_profiles(rbind(rise()), tt, "P")
  expect_true(all(is.na(e$dance_dlmo_run(pr, "fixed", list(threshold = 3))$phase_angle_h)))
  tb <- e$dance_dlmo_run(pr, "fixed", list(threshold = 3), bedtime = 23)
  expect_equal(tb$phase_angle_h, 2.5)
})

test_that("a dlmoR result is read back as hours on the profile's clock", {
  o <- as.POSIXct("2000-01-01", tz = "UTC")
  expect_equal(e$dance_dlmo_as_hours(o + 21.5 * 3600, o, c(18, 26)), 21.5)
  expect_equal(e$dance_dlmo_as_hours(o + 24.5 * 3600, o, c(18, 26)), 24.5)
  expect_equal(e$dance_dlmo_as_hours("00:30:00", o, c(18, 26)), 24.5)
  expect_equal(e$dance_dlmo_as_hours(as.difftime(21.25, units = "hours"), o, c(18, 26)), 21.25)
  # dlmoR 0.2.0 returns a time of day (hms): 07:00 in an 18:00-07:00 profile is
  # the next morning, not the evening before the window
  expect_equal(e$dance_dlmo_as_hours(structure(7 * 3600, class = c("hms", "difftime"), units = "secs"),
                                     o, c(18, 31)), 31)
  expect_equal(e$dance_dlmo_as_hours(structure(20.5 * 3600, class = c("hms", "difftime"), units = "secs"),
                                     o, c(18, 31)), 20.5)
  # with the hms package itself: as.numeric(units = "hours") on an hms gives seconds
  skip_if_not_installed("hms")
  expect_equal(e$dance_dlmo_as_hours(hms::hms(hours = 21, minutes = 3), o, c(18, 26)), 21.05)
  expect_equal(e$dance_dlmo_as_hours(hms::hms(minutes = 42), o, c(10, 32)), 24.7)
  # a date-time on some other day is a time of day
  expect_equal(e$dance_dlmo_as_hours(as.POSIXct("2024-05-03 21:00", tz = "UTC"), o, c(18, 26)), 21)
})

test_that("without dlmoR the hockey-stick is unavailable, never approximated", {
  skip_if(requireNamespace("dlmoR", quietly = TRUE))
  p <- list(id = "a", night = "1", group = NA, row = 1, t = tt, y = rise())
  r <- e$dance_dlmo_estimate(p, "hockey", list(threshold = 3))
  expect_identical(r$status, "unavailable"); expect_true(is.na(r$dlmo_h))
})

test_that("hockey-stick through dlmoR finds an onset near a known one", {
  skip_if_not_installed("dlmoR")
  set.seed(3); t5 <- seq(17, 25, by = 0.5)
  y <- 1 + 0.2 * rnorm(length(t5)) + pmax(0, t5 - 20.5)^1.6 * 3
  r <- e$dance_dlmo_hockey(t5, y, 3)
  expect_identical(r$status, "determined")
  expect_true(abs(r$time - 20.5) < 1)
})

sim_nights <- function(n = 14, shift = 1, groups = FALSE, seed = 1) {
  set.seed(seed)
  on1 <- 20.5 + rnorm(n, 0, 0.6); on2 <- on1 + shift + rnorm(n, 0, 0.3)
  g <- rep(c("A", "B"), length.out = n)
  if (groups) on2 <- on2 + ifelse(g == "B", 1, 0)
  t5 <- seq(18, 26, by = 0.5)
  f <- function(on) 1 + 6 * pmax(0, t5 - on) + rnorm(length(t5), 0, 0.05)
  Y <- t(sapply(seq_len(n), function(i) c(f(on1[i]), f(on2[i]))))
  pr <- e$dance_dlmo_profiles(Y, c(t5, t5 + 24), sprintf("S%02d", 1:n), group = g)
  list(tab = e$dance_dlmo_run(pr, "fixed", list(threshold = 3)), on1 = on1, on2 = on2)
}

test_that("the onsets of simulated profiles are recovered", {
  s <- sim_nights()
  n1 <- s$tab[s$tab$night == "Night 1", ]
  expect_true(all(n1$status == "determined"))
  # threshold 3 is reached 1/3 h after the rise starts
  expect_lt(max(abs(n1$dlmo_h - (s$on1 + 1 / 3))), 0.05)
})

test_that("nights compared within participant: paired t, and the mixed model with groups", {
  s <- sim_nights(shift = 1)
  cp <- e$dance_dlmo_compare(s$tab, within = TRUE)
  expect_true(cp$ok); expect_identical(cp$method, "paired t test")
  expect_lt(abs(cp$tests$estimate - 1), 0.25)
  expect_equal(cp$n_pairs, 14)
  # a participant without an onset on one night drops from the pairs, and is counted
  tb <- s$tab; tb$status[tb$id == "S01" & tb$night == "Night 2"] <- "no_onset"
  cp2 <- e$dance_dlmo_compare(tb, within = TRUE)
  expect_equal(cp2$n_pairs, 13); expect_equal(cp2$n_unpaired, 1)
  ds <- e$dance_dlmo_describe(tb, "dlmo_h", "night")
  expect_equal(ds$n_after[ds$cell == "Night 2"], 1)
  skip_if_not_installed("lmerTest")
  sg <- sim_nights(n = 20, shift = 1, groups = TRUE, seed = 2)
  cg <- e$dance_dlmo_compare(sg$tab, between = TRUE, within = TRUE)
  expect_true(cg$ok)
  expect_true(all(c("group", "night", "group:night") %in% cg$tests$effect))
  ct <- cg$tests[grepl("in B$", cg$tests$effect), ]
  expect_lt(abs(ct$estimate - 2), 0.4)
  expect_lt(cg$tests$p[cg$tests$effect == "group:night"], 0.01)
  # groups alone with two profiles per participant is refused
  expect_false(e$dance_dlmo_compare(sg$tab, between = TRUE)$ok)
})

test_that("two methods on the same profiles: their difference", {
  pr <- e$dance_dlmo_profiles(rbind(c(1, 1.2, 0.8, 1, 2, 4, 8, 12, 14)), tt, "P")
  tab <- e$dance_dlmo_run(pr, c("fixed", "variable"), list(threshold = 3, n_base = 3, k_sd = 2))
  ag <- e$dance_dlmo_agreement(tab)
  expect_equal(ag$n_both, 1)
  expect_equal(ag$mean_diff_h, 21.4 - 22.5)
})

test_that("a rise that starts between samples is interpolated early (sparse sampling)", {
  # the rise starts at 20:30 and reaches 3 at 20:50; hourly samples see 1 at
  # 20:00 and 4 at 21:00, so the straight line between them crosses 3 at 20:40
  y <- 1 + 6 * pmax(0, tt - 20.5)
  r <- e$dance_dlmo_cross(tt, y, 3)
  expect_equal(r$time, 20 + 2 / 3)
  expect_lt(r$time, 20.5 + 1 / 3)
})

test_that("three conditions within participants: the mixed model and every pair, Holm-adjusted", {
  skip_if_not_installed("lmerTest")
  set.seed(4); n <- 16; t5 <- seq(18, 26, by = 0.5)
  on <- 20.5 + rnorm(n, 0, 0.6)
  shift <- c(base = 0, sd = 1.5, rec = 0.5)            # condition effects on the onset, h
  f <- function(o) 1 + 6 * pmax(0, t5 - o) + rnorm(length(t5), 0, 0.05)
  Y <- do.call(rbind, lapply(names(shift), function(k) t(sapply(on + shift[[k]] + rnorm(n, 0, 0.2), f))))
  pr <- e$dance_dlmo_profiles(Y, t5, rep(sprintf("S%02d", 1:n), 3), night = rep(names(shift), each = n))
  tab <- e$dance_dlmo_run(pr, "fixed", list(threshold = 3))
  cp <- e$dance_dlmo_compare(tab, within = TRUE)
  expect_true(cp$ok)
  expect_match(cp$method, "Holm")
  pw <- cp$tests[cp$tests$test == "contrast (Satterthwaite)", ]
  expect_equal(nrow(pw), 3L)
  lv <- sort(names(shift))                              # factor levels: base, rec, sd
  want <- c(shift[[lv[2]]] - shift[[lv[1]]], shift[[lv[3]]] - shift[[lv[1]]], shift[[lv[3]]] - shift[[lv[2]]])
  expect_equal(pw$effect, c(sprintf("%s - %s", lv[2], lv[1]), sprintf("%s - %s", lv[3], lv[1]), sprintf("%s - %s", lv[3], lv[2])))
  expect_lt(max(abs(pw$estimate - want)), 0.2)
  expect_equal(pw$p_adj, stats::p.adjust(pw$p, "holm"))
  expect_true(all(pw$n == n))                           # every participant in both conditions of each pair
  expect_equal(cp$n_complete, n)
  expect_match(e$dance_dlmo_design_text(tab, within = TRUE), "3 nights / conditions \\(base, rec, sd\\); 16 participants")
  # one identifier per row (no repeat): nothing to pair, and said so
  tab2 <- tab; tab2$id <- paste0(tab2$id, "_", tab2$night)
  r2 <- e$dance_dlmo_compare(tab2, within = TRUE)
  expect_false(r2$ok)
  expect_match(r2$message, "must repeat")
})
