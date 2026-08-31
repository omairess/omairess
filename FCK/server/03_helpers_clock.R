# ==============================================================================
# server/03_helpers_clock.R — genuinely parsed clock times
#
# WHY THIS EXISTS
# ---------------
# WaPaa has two different notions of "time" and it is easy to mistake one for
# the other:
#
#   extract_time_values()      returns 1:n. It does NOT read the column names —
#                              its own comment says "Use sequential numbering to
#                              maintain chronology". This is what feeds
#                              values$time_numeric and every plot's x coordinate.
#   extract_hour_from_colname()  the real parser: Base9h -> 9, Base7h30 -> 7.5,
#                              KSS_9u_dag1 -> 9, X8.25 -> 8.25. It is used only
#                              to make axis LABELS, via get_hour_labels().
#
# Both are left exactly as they are, so every ported plot keeps the coordinates
# it had in the source app. This file adds a third, explicitly-named thing:
# clock times that are only reported when they can actually be trusted, for the
# places that need real elapsed time rather than a column counter (the cosinor
# tab's shared-time option, and the opt-in real-time smoothing).
# ==============================================================================

# Per-column clock hours, or NULL when the column names do not give trustworthy
# ones. "Trustworthy" = every column parses to a finite hour in [0, 24).
# extract_hour_from_colname() falls back to "any number in the name", which for
# labels like T1...T100 yields 1..100 — the range check rejects exactly that.
fck_clock_hours <- function(time_labels) {
  if (is.null(time_labels) || !length(time_labels)) return(NULL)
  hours <- suppressWarnings(vapply(time_labels, function(nm)
    tryCatch(as.numeric(extract_hour_from_colname(nm)),
             error = function(e) NA_real_),
    numeric(1), USE.NAMES = FALSE))
  if (any(!is.finite(hours))) return(NULL)
  if (any(hours < 0 | hours >= 24)) return(NULL)
  hours
}

# Elapsed hours since the first column, unwrapping midnight crossings, e.g.
# 22, 23, 0, 1 -> 0, 1, 2, 3. Same progression logic as WaPaa's
# calculate_time_positions(), which normalises the result to 0-1 for plotting;
# this returns it in hours, which is what a basis needs as argvals.
# NULL when fck_clock_hours() cannot vouch for the times.
fck_cumulative_hours <- function(time_labels) {
  hours <- fck_clock_hours(time_labels)
  if (is.null(hours)) return(NULL)
  n <- length(hours)
  if (n < 2) return(NULL)
  cum <- numeric(n)
  for (i in 2:n) {
    step <- if (hours[i] >= hours[i - 1]) hours[i] - hours[i - 1]
            else (24 - hours[i - 1]) + hours[i]
    cum[i] <- cum[i - 1] + step
  }
  # Two identical timestamps in a row would give a zero-width step and a basis
  # with duplicate breakpoints; that is a broken column set, not a time axis.
  if (any(diff(cum) <= 0)) return(NULL)
  cum
}

# TRUE when the columns are not evenly spaced in real time — the case where
# smoothing against the column index actually distorts the curves.
fck_spacing_is_uneven <- function(time_labels, tol = 1e-6) {
  cum <- fck_cumulative_hours(time_labels)
  if (is.null(cum)) return(FALSE)
  d <- diff(cum)
  (max(d) - min(d)) > tol * max(d)
}
