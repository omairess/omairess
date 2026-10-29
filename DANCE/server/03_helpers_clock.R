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

# EXCEL TIMES AS COLUMN NAMES. A header cell formatted as a time (08:00) is
# stored by Excel as a fraction of a day (1/3), and readxl hands the name over
# as that number, "0.33333333333333331", with "...5" appended when the same
# time occurs twice. The hour parser then read 0.333 h. This turns such names
# back into the clock time they display. Only when it is clearly what happened:
# at least two names that are numbers in [0, 1) falling on a whole minute of the
# day, at least one of them carrying the full binary precision Excel writes (a
# typed decimal hour such as 0.25 or 0.5 does not, and on its own is left as
# it is). A clock time that recurs is numbered, "08:00 #2", so the names stay
# unique; the parser reads "08:00 #2" as 8.
# Returns the names, converted where that applies, and which ones were.
dance_dayfraction_names <- function(nms) {
  nms <- as.character(nms)
  base <- sub("[.][.][.][0-9]+$", "", nms)
  x <- suppressWarnings(as.numeric(base))
  mins <- x * 1440
  cand <- grepl("^[0-9.]+([eE][-+]?[0-9]+)?$", base) & is.finite(x) & x >= 0 & x < 1 &
    abs(mins - round(mins)) < 1e-6
  long <- cand & nchar(gsub("[eE].*$", "", gsub("^0[.]?", "", base))) >= 10
  if (sum(cand) < 2 || !any(long)) return(list(names = nms, converted = rep(FALSE, length(nms))))
  m <- round(mins[cand])
  lab <- sprintf("%02d:%02d", m %/% 60, m %% 60)
  k <- stats::ave(seq_along(lab), lab, FUN = seq_along)
  lab[k > 1] <- sprintf("%s #%d", lab[k > 1], k[k > 1])
  out <- nms; out[cand] <- lab
  list(names = out, converted = cand)
}

# Per-column clock hours, or NULL when the column names do not give trustworthy
# ones. "Trustworthy" = every column parses to a finite hour in [0, 24).
# extract_hour_from_colname() falls back to "any number in the name", which for
# labels like T1...T100 yields 1..100 — the range check rejects exactly that.
dance_clock_hours <- function(time_labels) {
  if (is.null(time_labels) || !length(time_labels)) return(NULL)
  hours <- suppressWarnings(vapply(time_labels, function(nm)
    tryCatch(as.numeric(extract_hour_from_colname(nm)),
             error = function(e) NA_real_),
    numeric(1), USE.NAMES = FALSE))
  if (any(!is.finite(hours))) return(NULL)
  # "24" is midnight in many protocols (..., 22, 24, 2, 4, ...): read as 0 when
  # an earlier hour follows it, i.e. the recording crosses midnight there. A
  # series that merely ends at 24 (T1 ... T24) is still not taken for clock time.
  at24 <- which(hours == 24)
  if (length(at24) && all(vapply(at24, function(i)
        i < length(hours) && any(hours[(i + 1L):length(hours)] < 24), logical(1))))
    hours[at24] <- 0
  if (any(hours < 0 | hours >= 24)) return(NULL)
  hours
}

# Elapsed hours since the first column, unwrapping midnight crossings, e.g.
# 22, 23, 0, 1 -> 0, 1, 2, 3. Same progression logic as WaPaa's
# calculate_time_positions(), which normalises the result to 0-1 for plotting;
# this returns it in hours, which is what a basis needs as argvals.
# NULL when dance_clock_hours() cannot vouch for the times.
dance_cumulative_hours <- function(time_labels) {
  hours <- dance_clock_hours(time_labels)
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

# CLOCK TIMES IN RECORDING ORDER, for the cosinor. The columns are taken to be
# in chronological order, so a time EARLIER than the one before it is a
# midnight crossing: it and every later time move 24 h on. The cosinor used to
# apply this only to hand-entered times -- and by one PERIOD rather than 24 h,
# which is the same thing only when the period is 24 h -- while the clock times
# parsed at import (hours in [0, 24)) went in unwrapped, so a recording from
# 08:00 to 15:00 the next day put its second day BEFORE its first.
# Returns the times and whether anything was moved.
dance_unwrap_clock <- function(times, day = 24) {
  t <- as.numeric(times)
  if (length(t) > 1L) for (i in 2:length(t))
    if (is.finite(t[i]) && is.finite(t[i - 1]) && t[i] < t[i - 1])
      t[i:length(t)] <- t[i:length(t)] + day
  list(times = t, wrapped = !isTRUE(all.equal(t, as.numeric(times))))
}

# TRUE when the columns are not evenly spaced in real time — the case where
# smoothing against the column index actually distorts the curves.
dance_spacing_is_uneven <- function(time_labels, tol = 1e-6) {
  cum <- dance_cumulative_hours(time_labels)
  if (is.null(cum)) return(FALSE)
  d <- diff(cum)
  (max(d) - min(d)) > tol * max(d)
}
