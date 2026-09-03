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

# Data-driven default bandwidth (Taylor 2008, the circular plug-in rule):
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
  if (!is.finite(r_bar) || r_bar < 1e-6) return(period / 6)   # no scale: 4 h
  kappa <- if (r_bar < 0.53) 2 * r_bar + r_bar^3 + 5 * r_bar^5 / 6
           else if (r_bar < 0.85) -0.4 + 1.39 * r_bar + 0.43 / (1 - r_bar)
           else 1 / (r_bar^3 - 4 * r_bar^2 + 3 * r_bar)
  if (!is.finite(kappa) || kappa <= 0) return(period / 6)
  num <- 3 * n * kappa^2 * besselI(2 * kappa, nu = 2)
  den <- 4 * sqrt(pi) * besselI(kappa, nu = 1)^2
  kappa_bw <- if (is.finite(num) && is.finite(den) && den > 0) (num / den)^(2/5) else NA_real_
  if (!is.finite(kappa_bw) || kappa_bw <= 0) return(period / 6)
  bw <- (1 / sqrt(kappa_bw)) * period / (2 * pi)      # radians -> hours
  min(max(bw, period / 96), period / 4)               # 0.25 h .. 6 h on 24 h
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
