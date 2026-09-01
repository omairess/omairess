# ==============================================================================
# server/05_helpers_missing.R — where does every value in a smooth curve come from?
#
# The smoothing step fills gaps: it fits each subject's curve to the points that
# ARE there and then evaluates it at every time point, so a subject measured at
# 8 of 20 times comes out of the smoother with 20 values. Both source apps did
# this, correctly, and neither showed which 12 were invented.
#
# The distinction that matters is not observed-vs-missing but WHERE the missing
# value sits:
#
#   observed      a real measurement
#   interpolated  missing, but BETWEEN two observed points of that subject — the
#                 curve is supported on both sides, which is what smoothing is for
#   extrapolated  missing and OUTSIDE that subject's observed range — nothing
#                 constrains the curve there. A cubic spline evaluated far past
#                 its data does not level off; it follows its end polynomial.
#
# Both are "filled in", but only the first is interpolation in any defensible
# sense. Everything in this file exists to keep the two apart and visible.
# ==============================================================================

FCK_FILL_OBSERVED     <- 0L
FCK_FILL_INTERPOLATED <- 1L
FCK_FILL_EXTRAPOLATED <- 2L

FCK_FILL_LABELS <- c("Observed", "Interpolated", "Extrapolated")

# A status encoding, not a categorical one: "observed" is the quiet non-event and
# the two ways of being invented carry warning and critical. Validated for
# colour-vision deficiency (all-pairs deutan dE 18.0, normal-vision 18.7 against
# a #fcfcfb surface). Colour never carries this alone — the map has a legend and
# per-cell hover text, and the per-subject table gives the counts.
FCK_FILL_COLORS <- c(
  Observed     = "#e4e3dc",   # quiet neutral
  Interpolated = "#fab219",   # warning
  Extrapolated = "#d03b3b"    # critical
)

# Per-cell provenance for a subjects x time matrix of RAW (pre-smoothing) data.
# Returns an integer matrix of the three constants above.
fck_fill_status <- function(data_mat) {
  data_mat <- as.matrix(data_mat)
  n <- nrow(data_mat); p <- ncol(data_mat)
  obs <- !is.na(data_mat)
  status <- matrix(FCK_FILL_OBSERVED, n, p,
                   dimnames = list(rownames(data_mat), colnames(data_mat)))
  for (i in seq_len(n)) {
    which_obs <- which(obs[i, ])
    if (!length(which_obs)) {          # nothing at all: every value is invented
      status[i, ] <- FCK_FILL_EXTRAPOLATED
      next
    }
    inside <- seq_len(p) >= min(which_obs) & seq_len(p) <= max(which_obs)
    status[i, !obs[i, ] &  inside] <- FCK_FILL_INTERPOLATED
    status[i, !obs[i, ] & !inside] <- FCK_FILL_EXTRAPOLATED
  }
  status
}

# Per-subject counts, plus the two numbers that decide whether a curve is usable:
# how far the fit is carried past the data, and the widest gap it bridges.
fck_fill_per_subject <- function(status, hours_per_step = NA_real_) {
  n <- nrow(status)
  longest_run <- function(v) {
    if (!any(v)) return(0L)
    r <- rle(v); max(r$lengths[r$values])
  }
  out <- data.frame(
    Row          = seq_len(n),
    Observed     = rowSums(status == FCK_FILL_OBSERVED),
    Interpolated = rowSums(status == FCK_FILL_INTERPOLATED),
    Extrapolated = rowSums(status == FCK_FILL_EXTRAPOLATED),
    stringsAsFactors = FALSE
  )
  out$`% filled`     <- round(100 * (out$Interpolated + out$Extrapolated) / ncol(status), 1)
  out$`Longest gap`  <- vapply(seq_len(n), function(i)
    longest_run(status[i, ] == FCK_FILL_INTERPOLATED), integer(1))
  out$`First observed` <- vapply(seq_len(n), function(i) {
    w <- which(status[i, ] == FCK_FILL_OBSERVED); if (length(w)) w[1] else NA_integer_
  }, integer(1))
  out$`Last observed` <- vapply(seq_len(n), function(i) {
    w <- which(status[i, ] == FCK_FILL_OBSERVED); if (length(w)) w[length(w)] else NA_integer_
  }, integer(1))
  if (is.finite(hours_per_step)) {
    out$`Hours extrapolated` <- out$Extrapolated * hours_per_step
    out$`Longest gap (h)`    <- out$`Longest gap` * hours_per_step
  }
  out
}

# One line for the fit summary and the notification.
fck_fill_headline <- function(status) {
  tot <- length(status)
  sprintf("%d observed (%.0f%%), %d interpolated within the observed range (%.0f%%), %d extrapolated beyond it (%.0f%%)",
          sum(status == FCK_FILL_OBSERVED),     100 * mean(status == FCK_FILL_OBSERVED),
          sum(status == FCK_FILL_INTERPOLATED), 100 * mean(status == FCK_FILL_INTERPOLATED),
          sum(status == FCK_FILL_EXTRAPOLATED), 100 * mean(status == FCK_FILL_EXTRAPOLATED))
}

# Row labels for the map and the subject picker: the first scalar variable if
# there is one (a participant code is far more use than "row 47"), else the index.
fck_row_labels <- function(values) {
  n <- if (!is.null(values$data)) nrow(values$data) else 0L
  if (!n) return(character(0))
  lab <- NULL
  if (!is.null(values$covariates) && nrow(values$covariates) == n) {
    id_col <- names(values$covariates)[1]
    lab <- as.character(values$covariates[[id_col]])
    # repeated measures: a bare participant code repeats, so number the repeats
    if (anyDuplicated(lab)) lab <- paste0(lab, " #", ave(seq_along(lab), lab, FUN = seq_along))
  }
  if (is.null(lab) || length(lab) != n) lab <- paste("Row", seq_len(n))
  lab
}

# Mean spacing between time points in hours, when clock times could be parsed —
# turns "12 points extrapolated" into "24 hours extrapolated", which is the form
# a reader can judge.
fck_hours_per_step <- function(values) {
  cum <- tryCatch(fck_cumulative_hours(values$time_labels), error = function(e) NULL)
  if (is.null(cum) || length(cum) < 2) return(NA_real_)
  mean(diff(cum))
}
