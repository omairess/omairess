if (!exists("%||%", mode = "function")) `%||%` <- function(a, b) if (is.null(a)) b else a

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
# So there is now ONE rule, here. "No smoothing" means the observations joined
# by STRAIGHT LINES: an order-2 B-spline basis with a knot at every observation
# time, whose coefficients are the observed values themselves. The curve is the
# line plot of the data, and the mean curve is the line plot of the column
# means -- nothing is added between the measurement points.
#
# It used to be a cubic B-spline with nbasis = n_time on fda's default, EVENLY
# spaced knots. That also passed through every point, but its knots did not sit
# at the observation times, and between the points it swung far outside the
# data: on 18 hourly sleepiness ratings with means between 10 and 65 the mean
# curve ran from -2 to 79, with the largest swings next to the first and last
# observations.
#
# The penalised smoothing on the preprocessing tab is untouched: it remains
# the one place a roughness penalty is applied.
# ==============================================================================

# The straight-line basis: order 2 (piecewise linear) with a knot at every
# observation time, so there are exactly as many functions as observations and
# each one is 1 at its own time point and 0 at all the others.
dance_interpolating_basis <- function(argvals) {
  create.bspline.basis(rangeval = range(argvals), breaks = as.numeric(argvals), norder = 2)
}

# Build an fd object from a subjects x time matrix WITHOUT smoothing it: the
# coefficients ARE the observed values (no fitting step). argvals defaults to
# an even 0-1 grid, the range the fPCA family assumes. An fd object cannot
# carry NAs, so for this representation only a missing value is read off the
# straight line between the curve's neighbouring observations (held at the
# nearest observation before the first or after the last one); it used to be
# the participant's own mean, which put a step at every gap.
dance_fd_no_smoothing <- function(data_mat, argvals = NULL) {
  data_mat <- as.matrix(data_mat)
  n_time <- ncol(data_mat)
  if (is.null(argvals)) argvals <- seq(0, 1, length.out = n_time)
  for (i in seq_len(nrow(data_mat))) {
    ok <- !is.na(data_mat[i, ])
    if (!any(ok)) return(NULL)
    if (!all(ok))
      data_mat[i, !ok] <- if (sum(ok) == 1L) data_mat[i, ok]
                          else stats::approx(argvals[ok], data_mat[i, ok], argvals[!ok], rule = 2)$y
  }
  cf <- t(data_mat)
  dimnames(cf) <- list(NULL, rownames(data_mat))
  fd(cf, dance_interpolating_basis(argvals))
}

# Curves evaluated on a grid (time x curves), put back into an fd object on
# `basis` -- registered curves after time warping. On the straight-line basis
# of an unsmoothed run they are joined by straight lines through the grid
# values (no fitting, as above); on a smoothing basis they are projected onto
# it, as before. Every warped-curve path goes through here, so a run without a
# smoothing step is never smoothed after warping.
dance_is_straight_line_basis <- function(basis)
  identical(basis$type, "bspline") && basis$nbasis - length(basis$params) == 2L

dance_refd <- function(curves, time_points, basis = NULL) {
  curves <- as.matrix(curves)
  if (is.null(basis) || dance_is_straight_line_basis(basis))
    return(dance_fd_no_smoothing(t(curves), time_points))
  smooth.basis(time_points, curves, basis)$fd
}

# Used by the analysis tabs when values$fd_obj does not exist because the user
# went straight to an analysis. Builds the interpolating representation and
# SAYS SO, rather than quietly picking a smoothing basis on their behalf.
# Returns TRUE if an fd object is available afterwards.
dance_ensure_fd_obj <- function(values) {
  if (!is.null(values$fd_obj)) return(TRUE)
  dat <- if (!is.null(values$smooth_data)) values$smooth_data else values$data
  if (is.null(dat)) return(FALSE)

  fd <- tryCatch(dance_fd_no_smoothing(dat), error = function(e) NULL)
  if (is.null(fd)) {
    showNotification(
      "Could not build a curve representation from these data (a subject may be entirely missing). Apply smoothing on the Data Preprocessing tab first.",
      type = "error", duration = 12)
    return(FALSE)
  }
  values$fd_obj <- fd

  showNotification(
    paste("No smoothing step was applied, so each curve is your data points joined",
          "by straight lines: nothing has been smoothed. Visit 'Data",
          "Preprocessing/Smoothing' if you want a smooth."),
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
dance_auto_lambda <- function(data_mat, argvals, basisobj,
                            min_points_needed = 4,
                            log_range = c(-8, 4), n_grid = 25) {
  data_mat <- as.matrix(data_mat)
  rows <- which(rowSums(!is.na(data_mat)) >= min_points_needed)
  if (!length(rows)) return(NULL)

  # AUDIT (P12.1). This used to cap the objective at 60 subjects, taken as a
  # systematic subsample, under the argument that GCV is a smooth function of
  # the sample. The argument is not wrong, but it made the app's own description
  # of what it had done FALSE: the publication report said the parameter was
  # "chosen by generalised cross-validation over the whole sample", and on any
  # study with more than 60 usable subjects it was not. A method sentence that
  # overstates what was computed is the exact failure this whole audit exists to
  # prevent, and it is worse in the report than anywhere else, because that
  # sentence is the one that goes into a paper.
  #
  # So the cap is gone and the mean GCV is taken over every scorable subject.
  # Measured on 24 hourly columns with a 12-function B-spline basis, the search
  # is dominated by the fixed 25-point grid rather than by the sample:
  #     50 subjects (uncapped before and after)   4.0 s ->  4.2 s
  #    200 subjects   capped  5.3 s -> uncapped  17.3 s
  #    500 subjects   capped  4.4 s -> uncapped  42.0 s
  #   1000 subjects   capped  4.6 s -> uncapped  78.1 s
  #
  # Worth saying plainly: the cap was not producing a WRONG lambda. On the same
  # data the capped and uncapped searches returned 10.2 against 9.86 (200
  # subjects) and 11.4 against 9.84 (1000) -- differences far inside the range
  # over which GCV is flat. What it produced was a true number under a false
  # description. The caller reports n_used, so the cost is visible and the
  # sentence the report prints is now the sentence that is true.

  # AUDIT (P16). This called fda::smooth.basis ONCE PER SUBJECT, inside a loop
  # over a 25-point lambda grid: 25 x n fits to score one search. smooth.basis
  # accepts a MATRIX of curves and returns one GCV per column, and the only
  # thing that varies between subjects is which time points they have -- the
  # penalty, the basis and the argvals are shared. So subjects with the SAME
  # missingness pattern can be fitted in a single call.
  #
  # Complete data, the usual case, is one pattern and therefore one call.
  # Measured on 200 subjects x 24 time points with a 12-function B-spline:
  # 0.502 s per-subject against 0.003 s for the matrix call, 147x, with
  # max |difference| = 0 -- the same numbers, not an approximation.
  #
  # This is what pays for P12.1. Removing the 60-subject cap there was right --
  # the report claimed the whole sample and did not have it -- but it made every
  # automatic smoothing and every step of the n-basis sweep proportional to the
  # sample, and the user felt it. The cap stays gone; the redundant work goes
  # instead.
  pat <- apply(data_mat[rows, , drop = FALSE], 1, function(r) paste0(as.integer(!is.na(r)), collapse = ""))
  pat_groups <- split(rows, pat)

  mean_gcv <- function(log_lambda) {
    fdp <- fda::fdPar(basisobj, 2, 10^log_lambda)
    v <- unlist(lapply(pat_groups, function(idx) {
      ok <- !is.na(data_mat[idx[1], ])
      if (sum(ok) < 2) return(rep(NA_real_, length(idx)))
      ymat <- t(data_mat[idx, ok, drop = FALSE])       # time x subjects
      f <- tryCatch(fda::smooth.basis(argvals[ok], ymat, fdp),
                    error = function(e) NULL)
      if (is.null(f)) return(rep(NA_real_, length(idx)))
      g <- suppressWarnings(as.numeric(f$gcv))
      # one gcv per curve; anything else means the call did not do what we asked
      if (length(g) != length(idx)) return(rep(NA_real_, length(idx)))
      g[!is.finite(g)] <- NA_real_
      g
    }), use.names = FALSE)
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
# dance_n_harmonics(pca_res) -- how many components a pca.fd result actually holds
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
dance_n_harmonics <- function(pca_res) {
  if (is.null(pca_res) || is.null(pca_res$harmonics)) return(0L)
  h <- pca_res$harmonics
  n <- if (!is.null(h$coefs)) ncol(as.matrix(h$coefs)) else NA_integer_
  if (!is.finite(n) || n < 1) return(0L)
  as.integer(n)
}

# ==============================================================================
# dance_fpca(fd_obj, nharm, rotation) -- functional PCA, optionally rotated
# ==============================================================================
# fda::pca.fd(), then, on request, fda::varmx.pca.fd(): the VARIMAX rotation of
# the retained components (Ramsay & Silverman 2005, section 8.5), not
# reimplemented here. The rotation turns the nharm components so each has
# large values on fewer parts of the time axis, which can make them easier to
# name. What it changes, and the summary and report say so:
#   * the components are no longer ordered by variance, nor nested -- a rotated
#     set depends on how many components were retained;
#   * each component's share of the variance is recomputed (varprop), and the
#     total explained by the set is unchanged;
#   * the scores are rotated with them, so analyses of the scores use the
#     rotated components.
# The unrotated shares are kept (varprop_unrotated) so they can be reported.
dance_fpca <- function(fd_obj, nharm, rotation = c("none", "varimax")) {
  rotation <- match.arg(rotation)
  pr <- fda::pca.fd(fd_obj, nharm = nharm)
  pr$rotation <- "none"
  if (identical(rotation, "varimax") && nharm >= 2L) {
    vp0 <- pr$varprop
    pr <- fda::varmx.pca.fd(pr)
    pr$varprop_unrotated <- vp0
    pr$rotation <- "varimax"
  }
  pr
}

# ==============================================================================
# dance_smoothing_axis(input, values) / dance_smoothing_basis(...)
# ==============================================================================
# AUDIT (P9.3). The smoothing basis was built inline in server/20_smoothing.R,
# and the cross-validation diagnostic built its OWN with
#     create.bspline.basis(rangeval = c(1, n_time), nbasis = nb)
# on the integer column index. That is the production basis only in the
# default case. When the user has ticked cyclic smoothing, production fits a
# FOURIER basis; when they have ticked real clock times, production fits over
# elapsed hours with an uneven grid and a per-hour roughness penalty. In either
# case the CV's lambda is a penalty weight on a different penalty over a
# different basis, and the diagnostic's advice to type its smoothing factor
# into Data Preprocessing was wrong -- the same category of error as the
# mgcv-vs-fda ratio removed at P8.3, one level less obvious.
#
# Rather than adding a disclaimer, the two now build the SAME object. This
# repository has been bitten repeatedly by one statistical definition existing
# in two places and drifting; the fix is to have one.
#
# dance_smoothing_axis() resolves the time coordinates exactly as production
# does, and reports whether real clock time was actually usable.
dance_smoothing_axis <- function(input, values)
  dance_smoothing_axis_from(values$time_labels, ncol(values$data),
                            use_real_time = isTRUE(input$use_real_time),
                            cyclic = isTRUE(input$is_cyclic))

# The same axis from literal settings, so the exported script builds it with
# the app's own code (it has no `input` or `values`).
dance_smoothing_axis_from <- function(time_labels, n_time, use_real_time = FALSE,
                                      cyclic = FALSE) {
  real_time <- NULL
  if (isTRUE(use_real_time)) {
    real_time <- tryCatch(dance_cumulative_hours(time_labels),
                          error = function(e) NULL)
    if (is.null(real_time) || length(real_time) != n_time) real_time <- NULL
  }
  using_real_time <- !is.null(real_time)
  t_full <- if (using_real_time) real_time else seq_len(n_time)
  cyc <- dance_cyclic_period(time_labels, n_time, using_real_time)
  list(n_time = n_time, t_full = t_full, t_rng = range(t_full),
       using_real_time = using_real_time,
       cyclic = isTRUE(cyclic),
       period = cyc$period, period_source = cyc$source)
}

# ==============================================================================
# dance_smooth_curves() -- THE smoothing step, as a pure function
# ==============================================================================
# The computation of server/20_smoothing.R, moved here so the app and the
# exported script run the same code (emit_kernel). The script used to retype a
# simpler smooth -- every curve at once on 0-1, no per-curve handling of
# missing values, no real-time axis, no Fourier basis, no holding of the ends,
# no bounds -- so a run with any of those was not reproducible from it.
#
#   data      curves in rows, the observations as imported (NAs allowed)
#   axis      dance_smoothing_axis() / dance_smoothing_axis_from()
#   basis     dance_smoothing_basis(axis, nb, method)
#   method    "auto" (lambda by GCV, dance_auto_lambda), "manual" or "none"
#   lambda    the roughness penalty for "manual"
#   allow_extrapolation  FALSE holds each curve flat before its first and after
#             its last observation
#   bounds    NULL, or c(min, max) to clamp the smoothed values
#   progress  optional function(i, n) called as curves are smoothed
#
# Each curve is smoothed on its own observed points with fdPar(basis, 2,
# lambda) and evaluated at every time point (which is what fills its gaps); a
# curve the smoother cannot fit, or with fewer points than the basis needs
# (3 Fourier, 4 B-spline), is interpolated linearly and counted as failed. The
# fd object the fPCA family reads is the smoothed matrix re-expressed on 0-1
# with lambda = 0, on a basis of the same size and, for a Fourier basis, the
# same cycle (P19.3/P19.4).
dance_smooth_curves <- function(data, axis, basis, method = c("auto", "manual", "none"),
                                lambda = NULL, allow_extrapolation = FALSE,
                                bounds = NULL, progress = NULL) {
  method <- match.arg(method)
  data <- as.matrix(data)
  n_subjects <- nrow(data); n_time <- ncol(data)
  cyclic <- isTRUE(axis$cyclic)
  t_full <- axis$t_full; t_rng <- axis$t_rng
  time_points_01 <- if (isTRUE(axis$using_real_time)) (t_full - t_rng[1]) / diff(t_rng)
                    else seq(0, 1, length.out = n_time)
  if (identical(method, "none"))
    return(list(method = "none", smooth_mat = data, time_points_01 = time_points_01,
                fd_obj = dance_fd_no_smoothing(data, time_points_01)))

  min_points_needed <- if (cyclic) 3 else 4
  lambda_n_used <- NA_integer_; auto_failed <- FALSE
  lam <- if (identical(method, "auto")) {
    al <- dance_auto_lambda(data, t_full, basis, min_points_needed = min_points_needed)
    if (is.null(al)) { auto_failed <- TRUE; 0 } else { lambda_n_used <- al$n_used; al$lambda }
  } else as.numeric(lambda)
  fdParobj <- fda::fdPar(basis, 2, lam)

  smooth_mat <- matrix(NA, nrow = n_subjects, ncol = n_time)
  n_failed <- 0; failed_subjects <- c()
  df_vec <- rep(NA_real_, n_subjects); gcv_vec <- rep(NA_real_, n_subjects)
  for (i in seq_len(n_subjects)) {
    y_i <- data[i, ]
    valid_idx <- !is.na(y_i)
    n_valid <- sum(valid_idx)
    if (n_valid >= min_points_needed) {
      tryCatch({
        sb_i <- fda::smooth.basis(t_full[valid_idx], y_i[valid_idx], fdParobj)
        smooth_mat[i, ] <- as.vector(fda::eval.fd(t_full, sb_i$fd))
        df_vec[i] <- sb_i$df
        gcv_vec[i] <- sb_i$gcv
      }, error = function(e) {
        smooth_mat[i, ] <<- stats::approx(t_full[valid_idx], y_i[valid_idx],
                                          xout = t_full, rule = 2)$y
        n_failed <<- n_failed + 1
        failed_subjects <<- c(failed_subjects, i)
      })
    } else if (n_valid >= 2) {
      smooth_mat[i, ] <- tryCatch(
        stats::approx(t_full[valid_idx], y_i[valid_idx], xout = t_full, rule = 2)$y,
        error = function(e) rep(mean(y_i[valid_idx]), n_time))
    } else if (n_valid == 1) {
      smooth_mat[i, ] <- rep(as.numeric(y_i[valid_idx]), n_time)
      n_failed <- n_failed + 1
      failed_subjects <- c(failed_subjects, i)
    } else {
      n_failed <- n_failed + 1
      failed_subjects <- c(failed_subjects, i)
    }
    if (is.function(progress)) progress(i, n_subjects)
  }

  # outside each curve's observed range: held flat at its first/last
  # observation unless extrapolation was allowed
  n_clamped <- 0L
  if (!isTRUE(allow_extrapolation)) {
    for (i in seq_len(n_subjects)) {
      obs_i <- which(!is.na(data[i, ]))
      if (!length(obs_i)) next
      lo <- min(obs_i); hi <- max(obs_i)
      if (lo > 1) { smooth_mat[i, seq_len(lo - 1)] <- smooth_mat[i, lo]; n_clamped <- n_clamped + lo - 1L }
      if (hi < n_time) { smooth_mat[i, (hi + 1):n_time] <- smooth_mat[i, hi]; n_clamped <- n_clamped + n_time - hi }
    }
  }
  if (length(bounds) == 2L && all(is.finite(bounds)) && bounds[1] < bounds[2])
    smooth_mat <- pmin(pmax(smooth_mat, bounds[1]), bounds[2])

  # fit against the ORIGINAL values, where they were observed
  r_squared_vec <- numeric(n_subjects); rmse_vec <- numeric(n_subjects)
  for (i in seq_len(n_subjects)) {
    valid_idx <- !is.na(data[i, ])
    if (sum(valid_idx) > 1) {
      o <- data[i, valid_idx]; f <- smooth_mat[i, valid_idx]
      rmse_vec[i] <- sqrt(mean((o - f)^2))
      ss_tot <- sum((o - mean(o))^2)
      r_squared_vec[i] <- if (ss_tot > 0) 1 - sum((o - f)^2) / ss_tot else NA
    } else { r_squared_vec[i] <- NA; rmse_vec[i] <- NA }
  }

  # the fd object on 0-1: the same curves re-expressed, not smoothed again
  nb_used <- basis$nbasis
  basis_01 <- if (cyclic)
    fda::create.fourier.basis(rangeval = c(0, 1), nbasis = nb_used,
                              period = basis$params[1] / diff(t_rng))
  else fda::create.bspline.basis(rangeval = c(0, 1), nbasis = nb_used)
  fd_obj <- fda::smooth.basis(time_points_01, t(smooth_mat), fda::fdPar(basis_01, 2, 0))$fd

  list(method = method, smooth_mat = smooth_mat, fd_obj = fd_obj,
       time_points_01 = time_points_01, lambda = lam, lambda_n_used = lambda_n_used,
       auto_failed = auto_failed, nb_used = nb_used,
       cyclic_period = if (cyclic) basis$params[1] else NA_real_,
       df_vec = df_vec, gcv_vec = gcv_vec, r_squared = r_squared_vec, rmse = rmse_vec,
       n_failed = n_failed, failed_subjects = failed_subjects, n_clamped = n_clamped)
}

# ==============================================================================
# dance_cyclic_period(time_labels, n_time, using_real_time)
# ==============================================================================
# P19.3. The length of ONE CYCLE in the units of the smoothing axis, for the
# Fourier basis of cyclic smoothing.
#
# The basis used to be created without a period on the index axis, and fda then
# takes the WIDTH OF THE RANGE: columns 1..n span n - 1, so 24 hourly columns
# (00:00 ... 23:00) were fitted with a 23-column period. That forces the curve
# at 23:00 to equal the curve at 00:00, one hour early, and it misfits an exactly
# 24-hour signal by a visible amount at both ends. A cycle sampled at n evenly
# spaced columns without repeating its first point is n columns long, not n - 1.
#
# The rule, in order:
#   real clock time in use      24 (hours), as before;
#   clock hours parse AND the   24 / step columns, so 48 hourly columns are two
#   columns are evenly spaced   cycles of 24 and 00:00 ... 24:00 (25 columns,
#                               the endpoint repeated) is still 24;
#   otherwise                   n_time columns: the columns are taken to be one
#                               cycle, the last one a step before it repeats.
# `source` says which applied, so the summary can state it.
dance_cyclic_period <- function(time_labels, n_time, using_real_time = FALSE) {
  if (isTRUE(using_real_time))
    return(list(period = 24, source = "24 h of real clock time"))
  cum <- tryCatch(dance_cumulative_hours(time_labels), error = function(e) NULL)
  if (!is.null(cum) && length(cum) == n_time && n_time >= 2) {
    d <- diff(cum)
    if (all(is.finite(d)) && (max(d) - min(d)) <= 1e-6 * max(d)) {
      per <- 24 / mean(d)
      if (is.finite(per) && per >= 2)
        return(list(period = per,
                    source = sprintf("24 h = %s columns, from the clock times in the column names",
                                     format(signif(per, 6)))))
    }
  }
  list(period = n_time,
       source = sprintf("the %d columns taken as one cycle", n_time))
}

# The basis production fits. `nb` is the requested number of B-spline
# functions; it is ignored for a Fourier basis, which follows CIRCAREG's
# min(n_time, 13) rule (fda rounds an even count up to the next odd number,
# which is the behaviour CIRCAREG shipped).
# AUDIT (P10.2). P9.3 extracted the basis CONSTRUCTOR so production and the
# cross-validation diagnostic would build the same object, and left the basis
# COUNT rule duplicated: production capped at n_time, the CV at n_time - 2. On
# 16 time points with n_basis = 20 that is 16 functions against 14, so the claim
# that the CV fits "exactly the production basis" was not universally true --
# it held for cyclic smoothing, where the Fourier rule ignores the count, and
# for any series long enough for the cap not to bind. Extracting half a rule and
# leaving the other half duplicated is the same mistake this audit has made four
# times over; the whole rule lives here now.
dance_smoothing_nbasis <- function(input, n_time, method = NULL) {
  method <- method %||% (input$smooth_method %||% "auto")
  if (identical(method, "none")) return(min(20L, n_time - 2L))
  nb <- if (identical(method, "manual")) input$n_basis_manual else input$n_basis
  nb <- suppressWarnings(as.integer(nb))
  if (!length(nb) || !is.finite(nb)) nb <- 20L
  max(4L, min(nb, n_time))
}

# AUDIT (P11.3). `nb_fourier` exists for ONE caller: the GCV-vs-n-basis sweep,
# which has to vary the Fourier count to draw its curve. Everything else must
# leave it NULL and get the production rule. It is a parameter rather than a
# second constructor in the diagnostics file because the sweep still has to
# match production on everything ELSE -- the rangeval, the argvals and, for
# real-time data, period = 24 -- and building a second basis by hand is how the
# diagnostic came to fit a 13-hour period to 24-hour data in the first place.
#
# Note also what varying it means: for cyclic smoothing production PINS the
# count at min(n_time, 13) and ignores the user's n_basis entirely, so the
# sweep is exploratory there, not a recommendation. The sweep says so on screen.
dance_smoothing_basis <- function(axis, nb, method = "manual", nb_fourier = NULL) {
  nb <- max(4L, min(as.integer(nb), axis$n_time))
  if (isTRUE(axis$cyclic)) {
    nb_used <- if (is.null(nb_fourier)) min(axis$n_time, 13) else as.integer(nb_fourier)
    # P19.3: the period is stated, never left to fda's default (the width of the
    # range, which is one column short of a cycle). An axis built by hand
    # without one gets the same rule dance_cyclic_period() applies.
    per <- axis$period
    if (is.null(per) || !is.finite(per) || per <= 0)
      per <- if (isTRUE(axis$using_real_time)) 24 else axis$n_time
    b <- fda::create.fourier.basis(rangeval = axis$t_rng, nbasis = nb_used, period = per)
  } else if (identical(method, "none")) {
    # the straight-line representation dance_fd_no_smoothing() uses
    b <- fda::create.bspline.basis(rangeval = axis$t_rng, breaks = axis$t_full,
                                   norder = 2)
  } else {
    b <- fda::create.bspline.basis(rangeval = axis$t_rng, nbasis = nb)
  }
  b
}

# A one-line description, so a diagnostic can say which model its lambda
# belongs to instead of leaving the reader to assume.
dance_basis_label <- function(axis, basisobj) {
  sprintf("%s basis, %d functions, over %s [%.4g, %.4g]",
          if (isTRUE(axis$cyclic)) "Fourier" else "B-spline",
          basisobj$nbasis,
          if (isTRUE(axis$using_real_time)) "elapsed hours" else "column index",
          axis$t_rng[1], axis$t_rng[2])
}


# ==============================================================================
# dance_warp_amplitude(warping_results) -- one definition, used everywhere
# ==============================================================================
# AUDIT (P9.1). "How much was this curve warped" was computed in three places
# with three bodies. P8.2 corrected two of them and MISSED the third: the
# per-subject "Warping Amplitude Scores" table in server/40_fpca.R still read
#
#     if (!is.null(warp_results$warp_functions)) {
#       warp_amplitude[i] <- sqrt(mean((warp_functions[, i] - time_points)^2))
#     } else if (!is.null(warp_results$shifts)) { ... }
#
# and every linear-shift result carries warp_functions, so the first branch
# always won and the table went on showing the wrapped-shift artefact -- a
# PERIODIC ZERO SHIFT reported as 0.1 of the domain. Its third branch,
# abs(alpha_values - 1), was wrong for a different reason: since P4.1 only the
# power family has its identity at alpha = 1.
#
# The lesson is the recurrence, not the arithmetic. A duplicated statistical
# definition drifts; this is the fourth time in this audit. There is now one.
#
# A translation's amplitude is its displacement -- the shortest CIRCULAR one
# when the registration wrapped. An endpoint-preserving warp is measured by its
# RMS distance from the identity, which is what that means there.
dance_warp_amplitude <- function(warping_results) {
  if (is.null(warping_results)) return(NULL)
  sh <- warping_results$shifts
  wf <- warping_results$warp_functions
  tp <- warping_results$time_points

  if (identical(warping_results$method, "linear_shift") && !is.null(sh)) {
    sh <- as.numeric(sh)
    if (identical(warping_results$boundary, "periodic wrap")) {
      span <- if (!is.null(tp)) diff(range(tp)) else 1
      return(abs(((sh + span / 2) %% span) - span / 2))
    }
    return(abs(sh))
  }
  if (!is.null(wf)) {
    if (is.null(tp) || length(tp) != nrow(wf))
      tp <- seq(0, 1, length.out = nrow(wf))
    return(apply(wf, 2, function(h) sqrt(mean((h - tp)^2, na.rm = TRUE))))
  }
  NULL
}
