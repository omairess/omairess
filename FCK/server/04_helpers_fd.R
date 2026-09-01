# ==============================================================================
# server/04_helpers_fd.R — one rule for turning curves into an fd object
#
# THE PROBLEM THIS SOLVES
# -----------------------
# The fPCA / warping / fANOVA / clustering family cannot work on a matrix; it
# needs an fda `fd` object, and building one always means projecting onto a
# basis. In the source apps that projection was written in four places with
# two different sizes, and in two of them it fired without the user asking:
#
#   * "Raw data (no smoothing)" still built fd_obj on min(20, n_time - 2)
#     basis functions. On 24 hourly columns that is a real smooth — the label
#     said no smoothing, the object said otherwise.
#   * Running fPCA or fANOVA without visiting the smoothing tab built fd_obj
#     the same way, silently, from inside an analysis tab.
#
# So there is now ONE rule, here. "No smoothing" means an INTERPOLATING basis:
# nbasis = n_time, which makes smooth.basis a square system that reproduces
# the observed values at the observed times. The only thing lost is what is
# unavoidable in any fd representation — behaviour strictly between the
# measurement points.
#
# The penalised smoothing on the preprocessing tab is untouched: it remains
# the one place a roughness penalty is applied.
# ==============================================================================

# A basis that interpolates n_time points on `rangeval` rather than smoothing
# them. Cubic where there is room for it, lower order for very short series.
fck_interpolating_basis <- function(n_time, rangeval = c(0, 1)) {
  norder <- max(2L, min(4L, as.integer(n_time)))
  create.bspline.basis(rangeval = rangeval,
                       nbasis   = max(as.integer(n_time), norder),
                       norder   = norder)
}

# Build an fd object from a subjects x time matrix WITHOUT smoothing it.
# argvals defaults to an even 0-1 grid, the range the fPCA family assumes.
# NAs are mean-imputed for this representation only — an fd object cannot
# carry them — which is what the source apps did too.
fck_fd_no_smoothing <- function(data_mat, argvals = NULL) {
  data_mat <- as.matrix(data_mat)
  n_time <- ncol(data_mat)
  if (is.null(argvals)) argvals <- seq(0, 1, length.out = n_time)
  for (i in seq_len(nrow(data_mat))) {
    na_idx <- is.na(data_mat[i, ])
    if (any(na_idx)) {
      if (all(na_idx)) return(NULL)
      data_mat[i, na_idx] <- mean(data_mat[i, !na_idx], na.rm = TRUE)
    }
  }
  basis <- fck_interpolating_basis(n_time, range(argvals))
  smooth.basis(argvals, t(data_mat), fdPar(basis, 2, 0))$fd
}

# Used by the analysis tabs when values$fd_obj does not exist because the user
# went straight to an analysis. Builds the interpolating representation and
# SAYS SO, rather than quietly picking a smoothing basis on their behalf.
# Returns TRUE if an fd object is available afterwards.
fck_ensure_fd_obj <- function(values) {
  if (!is.null(values$fd_obj)) return(TRUE)
  dat <- if (!is.null(values$smooth_data)) values$smooth_data else values$data
  if (is.null(dat)) return(FALSE)

  fd <- tryCatch(fck_fd_no_smoothing(dat), error = function(e) NULL)
  if (is.null(fd)) {
    showNotification(
      "Could not build a curve representation from these data (a subject may be entirely missing). Apply smoothing on the Data Preprocessing tab first.",
      type = "error", duration = 12)
    return(FALSE)
  }
  values$fd_obj <- fd

  showNotification(
    paste("No smoothing step was applied, so these curves are represented by an",
          "INTERPOLATING basis: they pass exactly through your data points, and",
          "nothing has been smoothed. Visit 'Data Preprocessing/Smoothing' if you",
          "want a smooth."),
    type = "warning", duration = 12)
  TRUE
}
