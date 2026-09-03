# ==============================================================================
# server/07_helpers_circular.R — circular kernel density on a clock face
#
# Pure functions, no Shiny, so tests/circular_density_test.R can drive them.
#
# ORIENTATION
# -----------
# The plot puts DAY on top and NIGHT at the bottom, i.e. a 24-hour clock face
# with noon at the top and midnight at the bottom:
#
#         12:00                 06:00 -> left, 18:00 -> right, so the upper
#     06:00   18:00             half spans 06:00-18:00 (day) and the lower
#         00:00                 half 18:00-06:00 (night).
#
# plotly's angular axis is 0 degrees at east, counter-clockwise positive, so:
#
#     theta(h) = (270 - 15 * h) mod 360
#
# checked at the quarters: 00:00 -> 270 (bottom), 06:00 -> 180 (left),
# 12:00 -> 90 (top), 18:00 -> 0 (right). Hours therefore advance CLOCKWISE,
# which is what a reader expects of a clock face.
#
# DENSITY
# -------
# A von Mises kernel density, the circular analogue of a Gaussian KDE:
#
#     f(t) = 1 / (n * 2*pi * I0(kappa)) * sum_i exp(kappa * cos(t - t_i))
#
# It wraps by construction, so mass near midnight is not split between the two
# ends of a histogram the way a linear KDE of clock hours would split it.
# ==============================================================================

# Clock hour (in a period, default 24) -> plotly angular degrees, day at top.
fck_hour_to_theta <- function(hours, period = 24) {
  (270 - (360 / period) * (hours %% period)) %% 360
}

# ...and back, for hover text and tick labels.
fck_theta_to_hour <- function(theta_deg, period = 24) {
  ((270 - theta_deg) / (360 / period)) %% period
}

# Bandwidth as a von Mises concentration. `bw_hours` is the kernel's circular
# standard deviation in hours, which is the unit a chronobiologist can reason
# about; kappa ~= 1 / sigma^2 in radians.
fck_bw_to_kappa <- function(bw_hours, period = 24) {
  sigma_rad <- bw_hours * 2 * pi / period
  if (!is.finite(sigma_rad) || sigma_rad <= 0) return(NA_real_)
  1 / sigma_rad^2
}

# Data-driven default bandwidth (Taylor 2008, the circular plug-in rule), CAPPED.
#
# The cap is not cosmetic. The rule estimates concentration from R-bar, which is
# a GLOBAL quantity: for a symmetric bimodal sample R-bar is near zero however
# tight each mode is, so the rule concludes "uniform, smooth maximally" and
# returns a bandwidth that flattens the plot into a featureless disc exactly
# when there is structure worth seeing. Measured on simulated shapes, tight
# bimodal data drove it to the ceiling (density max/min 1.1 — a circle).
#
# Two fixes were tried. Deriving the concentration from higher trigonometric
# moments (the p-th moment picks up p-fold structure) rescued the bimodal and
# trimodal cases but INVENTED lobes in genuinely uniform data (max/min 6.8),
# which is the worse failure: a density plot that manufactures structure is
# worse than one that misses a subtle mode, and the bandwidth slider is right
# there. Capping at period/12 behaves sensibly across every shape tested
# (bimodal max/min 15.1, uniform 1.7) and can never manufacture anything.
#
#   kappa_bw = ( 3 n kappa^2 I2(2 kappa) / (4 sqrt(pi) I1(kappa)^2) ) ^ (2/5)
#
# with kappa the MLE-ish concentration of the sample (Best & Fisher's piecewise
# approximation, the same one watson_williams_test() uses). Returned in HOURS so
# it can seed the slider. Falls back to a wide-but-usable value when the sample
# is tiny or so dispersed that no scale is identified.
fck_default_bandwidth <- function(hours, period = 24, weights = NULL) {
  a <- (hours %% period) * 2 * pi / period
  ok <- is.finite(a)
  a <- a[ok]
  n <- length(a)
  if (n < 3) return(period / 12)                      # 2 h on a 24 h clock
  w <- if (is.null(weights)) rep(1, n) else {
    w <- weights[ok]; w[!is.finite(w) | w < 0] <- 0
    if (sum(w) <= 0) rep(1, n) else w
  }
  w <- w / sum(w)
  r_bar <- sqrt(sum(w * cos(a))^2 + sum(w * sin(a))^2)
  if (!is.finite(r_bar) || r_bar < 1e-6) return(period / 12)  # no scale: the cap
  kappa <- if (r_bar < 0.53) 2 * r_bar + r_bar^3 + 5 * r_bar^5 / 6
           else if (r_bar < 0.85) -0.4 + 1.39 * r_bar + 0.43 / (1 - r_bar)
           else 1 / (r_bar^3 - 4 * r_bar^2 + 3 * r_bar)
  if (!is.finite(kappa) || kappa <= 0) return(period / 12)
  # besselI overflows to Inf for large kappa, i.e. for the MOST concentrated
  # samples, and the rule then fell through to its fallback -- the widest
  # bandwidth, for the tightest data, which is backwards. The exponentially
  # scaled forms fix it exactly rather than approximately: I2(2k) carries
  # exp(2k) and I1(k)^2 carries exp(2k), so the factors cancel in the ratio.
  num <- 3 * n * kappa^2 * besselI(2 * kappa, nu = 2, expon.scaled = TRUE)
  den <- 4 * sqrt(pi) * besselI(kappa, nu = 1, expon.scaled = TRUE)^2
  kappa_bw <- if (is.finite(num) && is.finite(den) && den > 0) (num / den)^(2/5) else NA_real_
  if (!is.finite(kappa_bw) || kappa_bw <= 0) return(period / 12)
  bw <- (1 / sqrt(kappa_bw)) * period / (2 * pi)      # radians -> hours
  min(max(bw, period / 96), period / 12)              # 0.25 h .. 2 h on 24 h
}

# How peaked the shape is: max/min of the density. Near 1 means the ring is
# essentially round — worth saying out loud, because a round ring reads as a
# broken plot when it is in fact the correct answer for near-uniform data.
fck_density_contrast <- function(dens) {
  if (is.null(dens) || !length(dens$density)) return(NA_real_)
  mn <- min(dens$density)
  if (!is.finite(mn) || mn <= 0) return(Inf)
  max(dens$density) / mn
}

# The density itself, evaluated on a grid of clock hours.
# `weights` (e.g. amplitudes) let a strong rhythm count for more than a weak
# one; NULL gives every subject equal weight, which is what a Watson-Williams
# style analysis assumes.
fck_circular_density <- function(hours, bw_hours, period = 24, n_grid = 361,
                                 weights = NULL) {
  a <- (hours %% period) * 2 * pi / period
  ok <- is.finite(a)
  a <- a[ok]
  if (!length(a)) return(NULL)
  w <- if (is.null(weights)) rep(1, length(a)) else {
    w <- weights[ok]; w[!is.finite(w) | w < 0] <- 0
    if (sum(w) <= 0) rep(1, length(a)) else w
  }
  w <- w / sum(w)

  kappa <- fck_bw_to_kappa(bw_hours, period)
  if (!is.finite(kappa)) return(NULL)

  grid_h <- seq(0, period, length.out = n_grid)       # closes the ring
  grid_a <- grid_h * 2 * pi / period

  # exp(kappa*cos(.)) overflows for very small bandwidths; factor out the peak.
  dens <- vapply(grid_a, function(t)
    sum(w * exp(kappa * (cos(t - a) - 1))), numeric(1))
  dens <- dens / (2 * pi * besselI(kappa, nu = 0, expon.scaled = TRUE))

  list(hours = grid_h, density = dens, kappa = kappa, n = length(a))
}

# Mean direction and concentration, weighted the same way as the density, so
# the arrow drawn on the plot is the mean of what is actually plotted.
fck_circular_summary <- function(hours, period = 24, weights = NULL) {
  a <- (hours %% period) * 2 * pi / period
  ok <- is.finite(a)
  a <- a[ok]
  if (!length(a)) return(NULL)
  w <- if (is.null(weights)) rep(1, length(a)) else {
    w <- weights[ok]; w[!is.finite(w) | w < 0] <- 0
    if (sum(w) <= 0) rep(1, length(a)) else w
  }
  w <- w / sum(w)
  x <- sum(w * cos(a)); y <- sum(w * sin(a))
  r_bar <- sqrt(x^2 + y^2)
  list(mean_hour = (atan2(y, x) %% (2 * pi)) * period / (2 * pi),
       r_bar = r_bar, n = length(a))
}

# The night sector as plotly needs it: night wraps midnight, so it is drawn as
# one or two arcs depending on whether dusk < dawn.
fck_night_arcs <- function(dusk = 18, dawn = 6, period = 24) {
  dusk <- dusk %% period; dawn <- dawn %% period
  if (isTRUE(all.equal(dusk, dawn))) return(list())
  if (dusk < dawn) list(c(dusk, dawn)) else list(c(dusk, period), c(0, dawn))
}

# ==============================================================================
# Wrapping a recording onto one clock face
#
# For the polar PROFILE plot (radius = the measured value at each clock time,
# closed into a ring) as opposed to the acrophase density above.
#
# A protocol longer than the period visits the same clock time more than once —
# your 38 h sleep-deprivation runs hit 06:00 twice — so the columns have to be
# split by day before anything is drawn, or the ring doubles back on itself and
# the fill turns into knots.
# ==============================================================================

# Split elapsed hours into clock time and day number.
fck_wrap_to_clock <- function(cum_hours, period = 24) {
  list(clock = cum_hours %% period,
       day   = as.integer(floor(cum_hours / period)) + 1L)
}

# Put one day's points in clock order and close the ring by repeating the first
# point at the end, so a filled polygon has no seam.
fck_close_ring <- function(hours, values) {
  ok <- is.finite(hours) & is.finite(values)
  hours <- hours[ok]; values <- values[ok]
  if (length(hours) < 2) return(NULL)
  o <- order(hours)
  list(hours = c(hours[o], hours[o][1]), values = c(values[o], values[o][1]))
}

# A ribbon (mean +/- something) as ONE closed polygon: out along the upper edge,
# back along the lower edge reversed. plotly fills scatterpolar with 'toself',
# which needs exactly that ordering.
fck_band_ring <- function(hours, lo, hi) {
  ok <- is.finite(hours) & is.finite(lo) & is.finite(hi)
  hours <- hours[ok]; lo <- lo[ok]; hi <- hi[ok]
  if (length(hours) < 2) return(NULL)
  o <- order(hours)
  hours <- hours[o]; lo <- lo[o]; hi <- hi[o]
  list(hours  = c(hours, hours[1], rev(hours), hours[length(hours)]),
       values = c(hi,    hi[1],    rev(lo),    lo[length(lo)]))
}

# ==============================================================================
# The rhythmic part of a fitted cosinor, for the polar version of the fit plot
#
# predict_from_coefs() lays coefficients out as
#     c(mesor, [trend coefs...], beta_cos_1, beta_sin_1, beta_cos_2, ...)
# and finds the harmonics at an offset that depends on trend_type. Calling it
# with trend_type = "none" to drop a trend would therefore read the TREND
# coefficients as the first harmonic — a silent, plausible-looking wrong curve.
# This keeps the real trend_type for indexing and simply never adds the trend
# term, which is the arithmetic that plot does minus one line.
#
# Why drop the trend at all: a linear, log or saturating trend is not periodic,
# so it has no place on a clock face — 08:00 on the first day and 08:00 on the
# second are the same angle but different values, and the ring would not close.
# The polar plot therefore shows the rhythm; the trend stays on the fit plot.
# ==============================================================================
fck_rhythm_from_coefs <- function(coefs, time_vec, period, n_harmonics,
                                  trend_type = "none") {
  if (is.logical(trend_type)) trend_type <- if (trend_type) "linear" else "none"
  n_trend <- switch(as.character(trend_type),
                    "none" = 0, "linear" = 1, "log" = 1, "exp_sat" = 2, 0)
  off <- 1 + n_trend
  pred <- rep(coefs[1], length(time_vec))          # MESOR
  for (h in seq_len(n_harmonics)) {
    omega <- 2 * pi * h / period
    bc <- coefs[off + 2 * h - 1]
    bs <- coefs[off + 2 * h]
    if (!is.finite(bc) || !is.finite(bs)) next
    pred <- pred + bc * cos(omega * time_vec) + bs * sin(omega * time_vec)
  }
  pred
}

# ==============================================================================
# Night as a gradient rather than a block
#
# A hard-edged grey half-circle says "night starts at 18:00 and stops at 06:00",
# which is not what night does and not what a reader should take from a
# background. A ramp that fades in at dusk, is deepest around solar midnight and
# fades out at dawn carries the same information without asserting an edge.
#
# Returns one wedge per step: centre angle in HOURS, width in hours, and the
# colour. The window is followed FORWARD from dusk to dawn, so a window that
# crosses midnight (23:00 -> 07:00, the usual case) needs no special handling.
# ==============================================================================
fck_night_gradient <- function(dusk = 23, dawn = 7, period = 24, n = 72,
                               ramp = c("#ffffff", "#cde2fb", "#9ec5f4", "#6da7ec")) {
  dusk <- dusk %% period; dawn <- dawn %% period
  len <- (dawn - dusk) %% period
  if (!is.finite(len) || len <= 0) return(NULL)          # no night to draw
  n <- max(8L, as.integer(n))

  u <- (seq_len(n) - 0.5) / n                            # wedge centres in [0,1]
  # deepest at the middle of the window, white at both ends
  depth <- 1 - abs(2 * u - 1)
  cols <- grDevices::colorRampPalette(ramp)(101)[round(depth * 100) + 1]

  data.frame(hour  = (dusk + u * len) %% period,
             width = len / n,
             depth = depth,
             color = cols,
             stringsAsFactors = FALSE)
}
