# ==============================================================================
# tests/circular_density_test.R — the acrophase density and its clock face
#
# Two things are easy to get silently wrong here:
#   the ORIENTATION — day must be the upper half and night the lower, hours
#   running clockwise. A sign slip mirrors the plot, and a mirrored clock still
#   looks like a plausible clock.
#   the WRAPPING — a density that splits mass at midnight instead of wrapping is
#   exactly the failure a circular KDE exists to avoid, and on a ring it is not
#   obvious by eye.
#
# Run with:   Rscript tests/circular_density_test.R      (from the FCK directory)
# Needs no packages.
# ==============================================================================

app_dir <- if (dir.exists("server")) "." else "FCK"
eval(parse(file.path(app_dir, "server/07_helpers_circular.R"), encoding = "UTF-8"))

failures <- 0
check <- function(label, cond, detail = "") {
  if (isTRUE(cond)) cat(sprintf("ok  : %s\n", label))
  else { failures <<- failures + 1; cat(sprintf("FAIL: %s  %s\n", label, detail)) }
}
near <- function(a, b, tol = 1e-8) all(abs(a - b) < tol)

# --- orientation: plotly degrees are 0 = east, counter-clockwise -------------
check("midnight at the bottom", near(fck_hour_to_theta(0),  270))
check("06:00 on the left",      near(fck_hour_to_theta(6),  180))
check("noon at the top",        near(fck_hour_to_theta(12),  90))
check("18:00 on the right",     near(fck_hour_to_theta(18),   0))

# the upper half (theta in 0..180) must be exactly the daytime hours 06-18
mid_day   <- fck_hour_to_theta(seq(6.5, 17.5, by = 0.5))
mid_night <- fck_hour_to_theta(c(seq(18.5, 23.5, by = 0.5), seq(0.5, 5.5, by = 0.5)))
check("all daytime hours sit in the upper half", all(mid_day > 0 & mid_day < 180),
      sprintf("range %.0f-%.0f", min(mid_day), max(mid_day)))
check("all night hours sit in the lower half", all(mid_night > 180 & mid_night < 360),
      sprintf("range %.0f-%.0f", min(mid_night), max(mid_night)))

# hours must advance clockwise, i.e. decreasing plotly degrees
check("hours run clockwise", near(fck_hour_to_theta(13), fck_hour_to_theta(12) - 15))
check("round trip hour -> theta -> hour",
      near(fck_theta_to_hour(fck_hour_to_theta(c(0, 3.25, 9, 17.5, 23.75))),
           c(0, 3.25, 9, 17.5, 23.75), 1e-9))

# --- bandwidth ---------------------------------------------------------------
# a 1 h kernel on a 24 h clock is 2*pi/24 rad, so kappa = 1/sigma^2
check("bandwidth -> kappa", near(fck_bw_to_kappa(1), (24 / (2 * pi))^2, 1e-9))
check("wider bandwidth is less concentrated", fck_bw_to_kappa(3) < fck_bw_to_kappa(1))
tight <- fck_default_bandwidth(rep(14, 40) + rnorm(40, 0, 0.2))
loose <- fck_default_bandwidth(seq(0, 23, length.out = 40))
check("auto bandwidth is narrower for concentrated data", tight < loose,
      sprintf("tight %.2f vs loose %.2f", tight, loose))
check("auto bandwidth stays in range", tight > 0 && loose <= 6)
check("auto bandwidth survives n < 3", is.finite(fck_default_bandwidth(c(4, 5))))

# --- the density wraps -------------------------------------------------------
# acrophases side by side across midnight: the peak must be AT midnight, not
# split into two lumps at either end of the clock
d <- fck_circular_density(c(23, 23.5, 0, 0.5, 1), bw_hours = 1)
peak <- d$hours[which.max(d$density)]
check("wraps across midnight (peak at 00:00, not split)",
      min(abs(peak - 0), abs(peak - 24)) < 0.5, sprintf("peak at %.2f h", peak))
check("density is continuous around the seam",
      near(d$density[1], d$density[length(d$density)], 1e-9))
check("density integrates to ~1 over the circle",
      abs(sum(head(d$density, -1)) * (2 * pi / (length(d$density) - 1)) - 1) < 1e-3,
      sprintf("got %.4f", sum(head(d$density, -1)) * (2 * pi / (length(d$density) - 1))))

# a very small bandwidth must not overflow to NaN (exp(kappa) would)
tiny <- fck_circular_density(c(12, 12.1), bw_hours = 0.25)
check("survives a very small bandwidth", all(is.finite(tiny$density)))

# --- weighting ---------------------------------------------------------------
# one strong rhythm at 18:00 against three weak ones at 06:00
hrs <- c(18, 6, 6, 6); amp <- c(30, 1, 1, 1)
unw <- fck_circular_summary(hrs)
wtd <- fck_circular_summary(hrs, weights = amp)
check("unweighted mean follows the majority", abs(unw$mean_hour - 6) < 2,
      sprintf("%.2f h", unw$mean_hour))
check("amplitude weighting follows the strong rhythm", abs(wtd$mean_hour - 18) < 2,
      sprintf("%.2f h", wtd$mean_hour))

# R-bar behaves
check("R-bar ~ 1 for identical acrophases",
      fck_circular_summary(rep(9, 10))$r_bar > 0.999)
check("R-bar ~ 0 for evenly spread acrophases",
      fck_circular_summary(seq(0, 23, by = 3))$r_bar < 1e-9)

# --- the night sector wraps too ---------------------------------------------
arcs <- fck_night_arcs(18, 6)
check("night 18->06 is drawn as two arcs", length(arcs) == 2 &&
      near(arcs[[1]], c(18, 24)) && near(arcs[[2]], c(0, 6)))
check("a night that does not cross midnight is one arc",
      length(fck_night_arcs(1, 5)) == 1)
check("an empty night is no arcs", length(fck_night_arcs(6, 6)) == 0)

# --- a non-24 h period still maps sensibly -----------------------------------
check("period 12: half-period at the bottom", near(fck_hour_to_theta(0, 12), 270))
check("period 12: quarter-period on the left", near(fck_hour_to_theta(3, 12), 180))

if (failures) { cat("\n", failures, " failure(s).\n", sep = ""); quit(status = 1) }
cat("\nCircular density tests passed.\n")
