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

# ==============================================================================
# AUTOMATIC LAMBDA (P0.2)
# ==============================================================================
# The smoothing tab offered "Automatic smoothing (REML)". It set lambda = 0 and
# relied on a comment claiming smooth.basis() optimises internally. It does not:
# it smooths with the lambda it is handed and reports that lambda's GCV score.
# lambda = 0 is the unpenalised fit, and with nbasis capped at n_time that fit
# INTERPOLATES -- zero residual df, data reproduced to machine precision.
#
# This does the search that the label always promised, with the same estimator
# that performs the smoothing. The objective is the mean GCV over subjects: one
# lambda is chosen for the whole sample because a per-subject lambda would make
# curves incomparable, which is the entire premise of the fPCA family
# downstream. A coarse log grid locates the basin, then optimise() refines it --
# GCV in lambda is not reliably unimodal, so a bare optimise() over a wide
# interval can settle in a local minimum.
#
# Returns NULL when no subject has enough observed points to score.
fck_auto_lambda <- function(data_mat, argvals, basisobj,
                            min_points_needed = 4,
                            log_range = c(-8, 4), n_grid = 25) {
  data_mat <- as.matrix(data_mat)
  rows <- which(rowSums(!is.na(data_mat)) >= min_points_needed)
  if (!length(rows)) return(NULL)

  # Cap the work: GCV is a smooth function of the sample, so a large study does
  # not need every subject in the objective.
  if (length(rows) > 60) rows <- rows[round(seq(1, length(rows), length.out = 60))]

  mean_gcv <- function(log_lambda) {
    fdp <- fda::fdPar(basisobj, 2, 10^log_lambda)
    v <- vapply(rows, function(i) {
      ok <- !is.na(data_mat[i, ])
      f <- tryCatch(fda::smooth.basis(argvals[ok], data_mat[i, ok], fdp),
                    error = function(e) NULL)
      if (is.null(f)) return(NA_real_)
      g <- suppressWarnings(as.numeric(f$gcv))
      if (!length(g) || !is.finite(g[1])) NA_real_ else g[1]
    }, numeric(1))
    if (all(is.na(v))) return(Inf)
    mean(v, na.rm = TRUE)
  }

  grid <- seq(log_range[1], log_range[2], length.out = n_grid)
  scores <- vapply(grid, mean_gcv, numeric(1))
  if (all(!is.finite(scores))) return(NULL)
  k <- which.min(scores)
  lo <- grid[max(1, k - 1)]; hi <- grid[min(length(grid), k + 1)]
  best <- if (lo < hi)
    tryCatch(stats::optimize(mean_gcv, interval = c(lo, hi), tol = 1e-3)$minimum,
             error = function(e) grid[k])
  else grid[k]

  list(lambda = 10^best, log10_lambda = best, n_used = length(rows),
       gcv = mean_gcv(best))
}

# ==============================================================================
# fck_n_harmonics(pca_res) -- how many components a pca.fd result actually holds
# ==============================================================================
# AUDIT (P6.3). Four places asked this question as `length(pca_res$harmonics)`.
# `harmonics` is an **fd object**, and an fd object is a list of three elements
# (coefs, basis, fdnames) -- so that expression returns 3 for every PCA ever
# run, whatever nharm was. Every consumer took min(..., 3) and silently capped
# itself at three components: the loadings plot drew three, the effect-of-scores
# plot drew three, and the "Components to show" slider's ceiling was set to
# three, so a user who asked for five got three with nothing saying why.
#
# The count lives in the coefficient matrix: one column per harmonic.
fck_n_harmonics <- function(pca_res) {
  if (is.null(pca_res) || is.null(pca_res$harmonics)) return(0L)
  h <- pca_res$harmonics
  n <- if (!is.null(h$coefs)) ncol(as.matrix(h$coefs)) else NA_integer_
  if (!is.finite(n) || n < 1) return(0L)
  as.integer(n)
}
