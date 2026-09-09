# ==============================================================================
# server/01c_helpers_norm.R — the functional L2 norm
# ==============================================================================
# AUDIT (P20/R9). This used to live at the bottom of
# server/09_helpers_pcanova.R, which is sourced AFTER the modules that call it.
# Inside the running app that is harmless -- every server file is sourced into
# one environment before any handler runs -- but two consumers do not get that
# environment:
#
#   the exported analysis script, which is assembled from deparse()d kernels
#   plus their dependency closure; and
#
#   the pure-helper tests, which source individual server/*.R files.
#
# server/07_helpers_mixed_perm.R papered over the gap with
#
#     if (exists("dance_l2_norm", mode = "function")) dance_l2_norm(v, time)
#     else sum(diff(time) * (head(v,-1) + tail(v,-1)) / 2)
#
# and those two branches are NOT the same statistic: the first is
# sqrt(integral v^2 dt), the second is integral v dt. The global permutation
# p-value therefore depended on which files happened to be loaded -- measured on
# one fixture at a fixed seed, p = .09 with the helper present and p = .20
# without. Moving the definition ahead of every caller removes the branch, so
# there is one global statistic and it is the same one everywhere.
#
# ==============================================================================
# FUNCTIONAL L2 NORM (P1.3)
# ==============================================================================
# The fANOVA modules computed their global statistic as sqrt(sum(v^2)) over the
# evaluation grid. That is a vector norm, not a functional one: its value scales
# with how densely the grid happens to be sampled, so the same data evaluated on
# 50 versus 100 points gives different numbers, and on an unevenly spaced grid it
# silently weights the dense regions more.
#
# The L2 norm of a function is sqrt(integral v(t)^2 dt). Trapezoidal weights are
# exact for the piecewise-linear interpolant the app already draws, and reduce to
# the old constant-weight form (up to the grid spacing) when the grid is even.
dance_l2_norm <- function(v, argvals = NULL) {
  v <- as.numeric(v)
  ok <- is.finite(v)
  if (!any(ok)) return(NA_real_)
  if (is.null(argvals)) argvals <- seq(0, 1, length.out = length(v))
  argvals <- as.numeric(argvals)
  if (length(argvals) != length(v)) return(NA_real_)
  v <- v[ok]; argvals <- argvals[ok]
  if (length(v) < 2) return(abs(v[1]))
  o <- order(argvals); v <- v[o]; argvals <- argvals[o]
  d <- diff(argvals)
  w <- c(d[1], (utils::head(d, -1) + utils::tail(d, -1)), d[length(d)]) / 2   # trapezoid
  sqrt(sum(w * v^2))
}

