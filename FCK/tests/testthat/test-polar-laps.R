# ==============================================================================
# tests/testthat/test-polar-laps.R
#
# User report: "the harmonic fit and the circular plot of the harmonic fit
# don't match again."
#
# They do match -- the values are identical at every point, and the first test
# here proves it to machine precision. What differs is what a recording longer
# than 24 h LOOKS like on a 24 h dial: it laps, the passes overlap, and the
# resulting shape resembles nothing in the 2-D plot. So the fix is not to the
# arithmetic but to the display: laps can be separated, and the note reports
# the computed agreement instead of claiming it.
#
# The second half covers the related request: when the recording spans more
# than a period, folding it onto one dial AVERAGES the columns that share a
# clock time. Under a homeostatic rise 08:00 on day 1 and day 2 differ by the
# whole trend, so the average is a level that occurred on neither day. The
# spiral mode keeps every point.
# ==============================================================================

`%||%` <- function(a, b) if (is.null(a)) b else a
app_dir <- if (dir.exists("server")) "." else if (dir.exists("../../server")) "../.." else "FCK"
source(file.path(app_dir, "server/08_helpers_cosinor.R"))
source(file.path(app_dir, "server/03_helpers_clock.R"))
source(file.path(app_dir, "server/07_helpers_circular.R"))

# fck_clock_hours() calls the ported WaPaa column-name parser, which lives in
# the 18k-line server files. A stand-in for "KSS_8u00"-style labels keeps this
# test free of them.
extract_hour_from_colname <- function(nm) {
  m <- regmatches(nm, regexec("([0-9]+)u([0-9]+)?", nm))[[1]]
  if (length(m) >= 2 && nzchar(m[2]))
    as.numeric(m[2]) + (if (length(m) >= 3 && nzchar(m[3])) as.numeric(m[3]) / 60 else 0)
  else NA_real_
}

# Lift fck_fit_rings and fck_profile_rings out of the shipped file so the tests
# exercise the real code rather than a copy.
#
# Everything goes into ONE environment, and the helpers those functions call go
# in with them. Evaluating the definitions into the caller's frame instead
# leaves their closure environment without the stub above, and
# fck_profile_rings then silently returns NULL because fck_clock_hours cannot
# find the parser -- which looks exactly like a broken function.
.fenv <- new.env(parent = globalenv())
for (f in c("server/08_helpers_cosinor.R", "server/03_helpers_clock.R",
            "server/07_helpers_circular.R"))
  sys.source(file.path(app_dir, f), envir = .fenv)
assign("%||%", `%||%`, envir = .fenv)
assign("extract_hour_from_colname", extract_hour_from_colname, envir = .fenv)

.load_fn <- function(name) {
  src <- paste(readLines(file.path(app_dir, "server/74_polar_density.R"), warn = FALSE),
               collapse = "\n")
  i <- regexpr(paste0(name, " <- function"), src)
  stopifnot(i > 0)
  b <- substring(src, i)
  ch <- strsplit(b, "")[[1]]; d <- 0; started <- FALSE; e <- NA
  for (k in seq_along(ch)) {
    if (ch[k] == "{") { d <- d + 1; started <- TRUE }
    if (ch[k] == "}") { d <- d - 1; if (started && d == 0) { e <- k; break } }
  }
  eval(parse(text = substring(b, 1, e)), envir = .fenv)
  get(name, envir = .fenv)
}
fck_fit_rings     <- .load_fn("fck_fit_rings")
fck_profile_rings <- .load_fn("fck_profile_rings")

# predict_from_coefs, as the 2-D Fitted Curves plot calls it
pfc <- function(coefs, tv, period, nh, trend, t_offset = 0) {
  pred <- rep(coefs[1], length(tv))
  ntr <- switch(trend, "none" = 0, "linear" = 1, "log" = 1, "exp_sat" = 2, 0)
  off <- 1 + ntr
  if (trend != "none") pred <- pred + switch(trend,
    "linear" = coefs[2] * tv,
    "log" = coefs[2] * log(tv - t_offset + 1),
    "exp_sat" = coefs[2] * (1 - exp(-(tv - t_offset) / coefs[3])))
  for (h in seq_len(nh)) {
    w <- 2 * pi * h / period
    pred <- pred + coefs[off + 2 * h - 1] * cos(w * tv) + coefs[off + 2 * h] * sin(w * tv)
  }
  pred
}

mk_mod <- function(span = 26, start = 20) {
  tv <- seq(start, start + span, length.out = 20)
  G <- list(`1` = c(15, 78, 11, 17,  -9, 4, 3),
            `2` = c(18, 70, 13, 20, -12, 5, 2))
  gf <- lapply(names(G), function(g) list(
    mean_coefs = G[[g]], mean_mesor = G[[g]][1], n = 100,
    mean_amplitudes = c(sqrt(G[[g]][4]^2 + G[[g]][5]^2), sqrt(G[[g]][6]^2 + G[[g]][7]^2)),
    sd_amplitudes = c(3, 1)))
  names(gf) <- names(G)
  list(period = 24, n_harmonics = 2, trend_type = "exp_sat", time_vec = tv,
       t_offset = min(tv), t_center = mean(tv), group_fits = gf, coefs = G)
}

# ============================================================ the match itself
test_that("the polar fit equals the 2-D fit at every point, per group", {
  mod <- mk_mod(span = 26)
  fr <- fck_fit_rings(list(density_fit_scope = "recording", density_lap_mode = "spiral",
                           harmonic_show_ci = FALSE),
                      list(harmonic_model = mod))
  expect_false(is.null(fr))
  for (g in names(mod$coefs)) {
    rr <- fr$rings[[g]]
    t_abs <- seq(min(mod$time_vec), max(mod$time_vec), length.out = length(rr$pred))
    ref <- pfc(mod$coefs[[g]], t_abs, 24, 2, "exp_sat", mod$t_offset)
    expect_equal(rr$pred, ref, tolerance = 1e-12)
  }
})

test_that("the peak and trough clock times agree with the 2-D fit", {
  mod <- mk_mod(span = 26)
  fr <- fck_fit_rings(list(density_fit_scope = "recording", density_lap_mode = "spiral"),
                      list(harmonic_model = mod))
  for (g in names(mod$coefs)) {
    rr <- fr$rings[[g]]
    tf <- seq(min(mod$time_vec), max(mod$time_vec), length.out = 500)
    y <- pfc(mod$coefs[[g]], tf, 24, 2, "exp_sat", mod$t_offset)
    expect_equal(rr$hours[which.max(rr$pred)], tf[which.max(y)] %% 24, tolerance = 0.1)
    expect_equal(rr$hours[which.min(rr$pred)], tf[which.min(y)] %% 24, tolerance = 0.1)
  }
})

test_that("the agreement is computed and reported, not asserted", {
  fr <- fck_fit_rings(list(density_fit_scope = "recording", density_lap_mode = "spiral"),
                      list(harmonic_model = mk_mod(span = 26)))
  expect_true(is.finite(fr$agreement))
  expect_lt(fr$agreement, 1e-9)     # the number the note prints
})

# ==================================================================== lapping
test_that("a recording longer than a period is flagged as lapping", {
  expect_true(fck_fit_rings(list(density_fit_scope = "recording"),
                            list(harmonic_model = mk_mod(span = 26)))$spans_period)
  expect_false(fck_fit_rings(list(density_fit_scope = "recording"),
                             list(harmonic_model = mk_mod(span = 22)))$spans_period)
})

test_that("per_lap splits into one trace per day and loses nothing", {
  mod <- mk_mod(span = 30)      # 30 h: two laps
  sp <- fck_fit_rings(list(density_fit_scope = "recording", density_lap_mode = "spiral"),
                      list(harmonic_model = mod))
  pl <- fck_fit_rings(list(density_fit_scope = "recording", density_lap_mode = "per_lap"),
                      list(harmonic_model = mod))
  expect_equal(length(sp$rings), 2)                 # two groups
  expect_gt(length(pl$rings), length(sp$rings))     # split further, by day
  expect_true(all(grepl("day", names(pl$rings))))
  # every point survives the split
  for (g in c("1", "2")) {
    parts <- pl$rings[names(pl$rings)[vapply(pl$rings, function(z) identical(z$parent, g), TRUE)]]
    n_split <- sum(vapply(parts, function(z) length(z$pred), 1L))
    expect_equal(n_split, length(sp$rings[[g]]$pred))
  }
})

test_that("a recording inside one period is not split", {
  mod <- mk_mod(span = 20)
  pl <- fck_fit_rings(list(density_fit_scope = "recording", density_lap_mode = "per_lap"),
                      list(harmonic_model = mod))
  expect_equal(sort(names(pl$rings)), c("1", "2"))   # untouched
})

test_that("a lapping arc must not be sorted onto one dial", {
  mod <- mk_mod(span = 30)
  fr <- fck_fit_rings(list(density_fit_scope = "recording", density_lap_mode = "spiral"),
                      list(harmonic_model = mod))
  h <- fr$rings[["1"]]$hours
  # clock hours are NOT monotone -- they wrap -- which is the whole point;
  # sorting them would fold day 2 onto day 1
  expect_false(!is.unsorted(h))
  expect_true(any(diff(h) < 0))
})

# ============================================== the profile: averaging vs spiral
mk_vals <- function(n_sub = 30, seed = 5) {
  set.seed(seed)
  # 08:00 -> 08:00 the NEXT day plus 6 h: clock 08 is visited twice
  hrs <- c(8:23, 0:23, 0:6)
  labs <- sprintf("KSS_%du00", hrs)
  cum <- cumsum(c(0, rep(1, length(hrs) - 1)))
  Y <- t(vapply(seq_len(n_sub), function(i)
    20 + 2.5 * cum + 10 * sin(2 * pi * cum / 24) + rnorm(length(cum), 0, 1),
    numeric(length(cum))))
  list(smooth_data = Y, time_labels = labs, time_clock = hrs %% 24,
       covariates = NULL)
}

test_that("folding onto one dial averages the repeated clock times", {
  v <- mk_vals()
  pr <- fck_profile_rings(list(density_profile_mode = "average"), v)
  expect_false(is.null(pr))
  z <- pr[[1]]
  expect_false(isTRUE(z$spiral))
  expect_equal(length(z$hours), 24)          # folded to one dial
  expect_true(any(z$n_folded > 1))           # and some hours were averaged
})

test_that("the spiral keeps every point and averages nothing across days", {
  v <- mk_vals()
  pr <- fck_profile_rings(list(density_profile_mode = "spiral"), v)
  z <- pr[[1]]
  expect_true(isTRUE(z$spiral))
  expect_equal(length(z$hours), ncol(v$smooth_data))   # nothing folded away
  expect_equal(length(z$mean), ncol(v$smooth_data))
  # the same clock hour appears more than once, at different values
  rep_h <- z$hours[duplicated(z$hours)][1]
  vals <- z$mean[z$hours == rep_h]
  expect_gt(length(vals), 1)
  expect_gt(diff(range(vals)), 5)    # they differ by the trend, as they should
})

test_that("averaging reports a level that occurred on neither day", {
  # the substantive reason the default changed
  v <- mk_vals()
  av <- fck_profile_rings(list(density_profile_mode = "average"), v)[[1]]
  sp <- fck_profile_rings(list(density_profile_mode = "spiral"), v)[[1]]
  h <- 8
  folded <- av$mean[av$hours == h]
  actual <- sp$mean[sp$hours == h]
  expect_equal(length(folded), 1)
  expect_gt(length(actual), 1)
  expect_true(folded > min(actual) && folded < max(actual))   # between, not equal
  expect_false(any(abs(actual - folded) < 1e-6))
})

test_that("per_day gives one ring per day", {
  v <- mk_vals()
  pr <- fck_profile_rings(list(density_profile_mode = "per_day"), v)
  expect_gt(length(pr), 1)
  expect_true(all(grepl("day", names(pr))))
  expect_true(all(vapply(pr, function(z) !isTRUE(z$spiral), logical(1))))
})

test_that("the old checkbox still selects the old behaviour", {
  v <- mk_vals()
  a <- fck_profile_rings(list(density_avg_days = TRUE), v)
  d <- fck_profile_rings(list(density_avg_days = FALSE), v)
  expect_equal(length(a), 1)
  expect_gt(length(d), 1)
})

test_that("a within-day recording behaves the same in every mode", {
  set.seed(9)
  hrs <- 8:20
  v <- list(smooth_data = t(vapply(1:10, function(i)
              20 + seq_along(hrs) + rnorm(length(hrs), 0, .5), numeric(length(hrs)))),
            time_labels = sprintf("KSS_%du00", hrs), time_clock = hrs)
  for (m in c("spiral", "average", "per_day")) {
    pr <- fck_profile_rings(list(density_profile_mode = m), v)
    expect_false(is.null(pr), info = m)
    expect_equal(length(pr[[1]]$hours), length(hrs), info = m)
  }
})
