# ==============================================================================
# server/05_helpers_warp.R -- the curve-registration kernels, in one place
# ==============================================================================
# AUDIT (P11.1 and P11.2). Two problems with one cause: these four functions
# used to live inside server/40_fpca.R, and one of them was nested inside a
# reactive observer.
#
# P11.2 -- THE LIVE DEFECT. fck_landmark_warp() was defined at depth 11, inside
# the `else if (warping_method == "landmark")` branch of the observer, while its
# only caller, landmark_alignment_simple(), was defined at the top level of the
# same file. R scopes functions LEXICALLY, so the caller's enclosing environment
# is the server environment and the observer's frame is not on its search path:
# the call could never resolve. Every landmark registration therefore died with
#     Error in landmark_alignment: could not find function "fck_landmark_warp"
# which the surrounding tryCatch turned into NULL, which the caller treated as
# "warping failed" and REPLACED WITH A LINEAR SHIFT. So a user who selected
# landmark registration silently received shift registration instead, with no
# on-screen indication: the notice went to the console. This was introduced by
# the P5.5/P5.6 landmark work -- the function was written correctly and put in
# the wrong place. Moving it to a file that is sourced at top level fixes it,
# and tests/registration_kernel_test.R calls the shipped landmark path so that
# a definition which cannot be reached fails a test rather than a user.
#
# P11.1 -- WHY A SEPARATE FILE, AND NOT JUST A MOVE. server/90_export.R wrote
# its own copy of the registration algorithms with add("..."), and that copy had
# not been updated since before P0.8. Measured, live estimator versus exported
# script, on a periodic profile with known displacements:
#
#     true delta    live s     exported s
#      0.05        -0.0505       0.0050
#      0.10        -0.1010       0.0090
#     -0.10        +0.1010      -0.0080
#
# The exported script recovered about a twelfth of the shift AND REVERSED ITS
# SIGN (it kept the 0.1 and 0.5 attenuation constants P0.8 removed, and had no
# periodic branch at all); the registered curves differed by up to 2.17 on
# curves of range +/-1.4. The parametric export was equally stale: on curves
# needing NO warping it pinned the quadratic family at the range boundary and
# deformed the time axis by 0.4999 of the domain, the exact defect P4.1 fixed,
# while the live kernel returned 0.0003. The `logistic` family the UI offers was
# absent from the exported switch, so it fell through to the identity.
#
# The generator already had the machinery to prevent this -- emit_kernel() in
# 90_export.R deparses the app's OWN function objects and resolves their helper
# closure with codetools::findGlobals -- and sections 7, 10 and 11 use it. The
# warping section did not, because these functions were not reachable as plain
# top-level objects. They are now, so the export emits the real thing and the
# two cannot drift apart again. This is the fifth time in this audit that a
# duplicated definition has drifted; the fix is always the same, and it is the
# reason this file exists rather than a comment telling the next person to keep
# the copies in step.
#
# CONTRACT (module-wide, unchanged from P5.6):
#   h maps REGISTERED time to ORIGINAL time, and registration is always
#       registered_i <- approx(time_points, original_i, xout = h_i)$y
#
# Everything here is PURE: no `input`, no `values`, no session, no notification.
# That is what makes it emittable, testable and shared. landmark_alignment_simple
# therefore takes landmark_points as an argument and RETURNS its rejection
# message in `rejection_warning` instead of calling showNotification(); the
# observer in 40_fpca.R shows it.
# ==============================================================================

# ==========================================================================
# fck_landmark_warp(ref, own, time_points)
#
# The one place a landmark warp is built. Returns h with the module-wide
# contract (P5.6):
#
#     h maps REGISTERED time to ORIGINAL time, and registration is always
#         registered_i <- approx(time_points, original_i, xout = h_i)$y
#
# so h(t) answers "which point of the original curve belongs at registered
# time t". Every method in this file returns h in that direction. Before P5.6
# the two landmark branches disagreed: the manual branch used this contract,
# while the automatic branch built h from (own -> ref) and then applied it as
# approx(h, curve, xout = t), which is the INVERSE map. The registered curves
# came out plausible either way, because the second form inverts the map
# while interpolating -- but warp_functions then held two different objects
# depending on which branch ran, and every downstream statistic that reads
# warp_functions - time_points (the phase metrics, the warping plots, the
# per-subject table) was comparing incomparable things across methods.
#
# AUDIT (P5.5): the manual branch also had no monotonicity requirement at
# all. It took whatever the peak/valley search returned, pasted 0 and 1 on
# the ends and interpolated. On a noisy curve the detected landmarks can
# cross -- the search alternates peaks and valleys in fixed windows, and
# nothing stopped landmark 3 landing before landmark 2 -- and a crossed pair
# makes h non-monotone, which is not a reparameterisation of time. It is now
# a hard precondition: crossed or duplicated landmarks are REJECTED and the
# subject falls back to the identity, reported, rather than being registered
# with a fold in it.
#
# Returns NULL when the landmarks do not define a valid warp.
# ==========================================================================
fck_landmark_warp <- function(ref, own, time_points) {
  if (length(ref) != length(own) || !length(ref)) return(NULL)
  if (any(!is.finite(ref)) || any(!is.finite(own))) return(NULL)

  t0 <- time_points[1]; t1 <- time_points[length(time_points)]
  # order both by the REFERENCE position, so a reference given out of order
  # does not silently reorder the correspondence
  o   <- order(ref)
  ref <- ref[o]; own <- own[o]

  # strictly increasing, strictly interior, on BOTH sides
  tol <- 1e-6 * max(1, t1 - t0)
  if (any(diff(ref) <= tol) || any(diff(own) <= tol)) return(NULL)
  if (ref[1] <= t0 + tol || ref[length(ref)] >= t1 - tol) return(NULL)
  if (own[1] <= t0 + tol || own[length(own)] >= t1 - tol) return(NULL)

  kx <- c(t0, ref, t1)     # registered time
  ky <- c(t0, own, t1)     # original time
  h  <- stats::approx(kx, ky, xout = time_points, rule = 2)$y

  # the invariants this function exists to guarantee
  if (!isTRUE(all.equal(h[1], t0, tolerance = 1e-8))) return(NULL)
  if (!isTRUE(all.equal(h[length(h)], t1, tolerance = 1e-8))) return(NULL)
  if (any(diff(h) <= 0)) return(NULL)
  h
}


# Warping functions with better error handling
linear_shift_alignment <- function(fd_obj, periodic = FALSE, 
                                   allow_dilation = FALSE, 
                                   dilation_range = c(0.95, 1.05),
                                   reference = "mean",
                                   time_points = seq(0, 1, length.out = 100)) {
  
  tryCatch({
    n_curves <- ncol(fd_obj$coefs)
    n_time <- length(time_points)
    
    # Evaluate curves
    curves <- eval.fd(time_points, fd_obj)
    
    # Validate curves
    if(is.null(curves) || ncol(curves) == 0) {
      stop("No valid curves to align")
    }
    
    # Get reference curve
    if(reference == "mean") {
      ref_curve <- rowMeans(curves)
    } else if(reference == "median") {
      ref_curve <- apply(curves, 1, median)
    } else {
      ref_curve <- curves[,1]
    }
    
    # Initialize warping
    warp_functions <- matrix(NA, n_time, n_curves)
    registered_curves <- matrix(NA, n_time, n_curves)
    shifts <- numeric(n_curves)
    extrap_frac <- numeric(n_curves)   # P4.2: how much of each warp left the domain

    # Perform alignment for each curve
    for(i in 1:n_curves) {
      # Find best shift using cross-correlation
      if(periodic) {
        # AUDIT (P10.1): the circular estimator was off by exactly one FFT bin,
        # and it was fed a duplicated endpoint. Two separate defects.
        #
        #   1. R's which.max() is 1-BASED, and element 1 of the inverse FFT is
        #      ZERO lag. The old code took `max_idx` as the lag directly, so
        #      every estimate was displaced by one grid step. Measured on the
        #      default 100-point grid with two IDENTICAL curves: s = -0.0100,
        #      i.e. a curve needing no shift was moved by 1% of the cycle --
        #      14.4 minutes on a 24-hour day. Every other estimate carried the
        #      same offset (a true 0.05 delay was estimated as 0.06).
        #
        #   2. The evaluation grid is seq(0, 1, length.out = n_time), which
        #      contains BOTH endpoints. On a cycle those are the same phase, so
        #      a circular cross-correlation must not see both: the duplicated
        #      sample biases the correlation and makes the lag resolution
        #      1/n_time where it should be 1/(n_time - 1).
        #
        # This is upstream of everything the registration produces. The shift
        # is APPLIED to the curve, so the per-subject R-squared, RMSE and
        # correlation in the warping table were computed on a curve that had
        # been moved when it should not have been: on the zero-shift case the
        # old code reported R2 = .9955, RMSE = .0557, r = .9978 for a curve
        # identical to the reference. With the fix those are exactly 1, 0, 1.
        #
        # Verified across true displacements 0, +/-0.05, +/-0.10 and 0.25: the
        # corrected estimator recovers each to within 1/(n_time - 1), which is
        # the grid resolution and not a systematic error.
        circ_idx <- seq_len(n_time - 1L)      # drop the duplicated endpoint
        n_circ   <- length(circ_idx)
        curve_c  <- curves[circ_idx, i]
        ref_c    <- ref_curve[circ_idx]

        fft_curve  <- fft(curve_c - mean(curve_c))
        fft_ref    <- fft(ref_c - mean(ref_c))
        cross_corr <- Re(fft(Conj(fft_ref) * fft_curve, inverse = TRUE)) / n_circ

        lag_idx <- which.max(cross_corr) - 1L   # bin 1 is zero lag
        if (lag_idx > n_circ / 2) lag_idx <- lag_idx - n_circ
        shifts[i] <- -lag_idx / n_circ
      } else {
        # Standard cross-correlation
        ccf_result <- ccf(curves[,i], ref_curve, lag.max = floor(n_time/4), 
                          plot = FALSE, na.action = na.pass)
        if(!is.null(ccf_result$acf) && length(ccf_result$acf) > 0) {
          best_lag <- ccf_result$lag[which.max(ccf_result$acf)]
          # AUDIT (P0.8): was  best_lag / n_time * 0.1  -- the estimated lag
          # was multiplied by 0.1 here and by a further 0.5 below, so the warp
          # carried 5% of the shift it had just measured. Neither constant was
          # justified anywhere. The measured lag is now used as measured.
          shifts[i] <- best_lag / n_time
        } else {
          shifts[i] <- 0
        }
      }
      
      # ---- AUDIT (P0.8): the warp is the shift, and nothing else -----------
      # This block used to read:
      #     base_warp  <- time_points - shifts[i] * 0.5
      #     distortion <- sin(pi * time_points) * runif(1, -0.03, 0.03)
      #     warp_functions[,i] <- pmin(1, pmax(0, base_warp + distortion))
      # under the comment "# Add slight S-curve for visualization". A random
      # draw was being added to an ESTIMATED transformation, on a module with
      # no set.seed, so the same data gave different answers on every run --
      # measured at up to 5.9x the size of the shift being estimated. The
      # clipping then destroyed monotonicity at the ends: a clipped warp maps
      # a whole interval of new time onto one old time, and approx() hands
      # back the same value for all of it.
      #
      # A shift warp is h(t) = t - s. It is monotone by construction and
      # deterministic. Forcing the endpoints by ASSIGNMENT (as the pre-P0.8
      # code did) is what broke monotonicity when |s| was large, so the shift
      # is limited to what the domain can absorb instead. See P4.2 below for
      # what h actually maps onto, and for the boundary rule.
      # AUDIT (P4.2): the comment above used to claim three properties for
      # this warp -- monotone by construction, deterministic, and (in its own
      # words) anchored at the endpoints. The first two are true. THE THIRD IS
      # FALSE, and a reviewer was right to call it out. (The exact old phrase
      # is deliberately not repeated here: tests/warp_family_test.R greps the
      # source for it, and a fix that quotes the sentence it is removing makes
      # its own guard vacuous -- which has happened twice in this audit.)
      # h(t) = t - s maps [0, 1] onto [-s, 1 - s]. For s = 0.1 that is
      # h(0) = -0.1 and h(1) = 0.9: it is not a map of [0,1] to itself and it
      # does not fix the endpoints. It is a TRANSLATION, which is what shift
      # registration is (Ramsay & Silverman's shift model), and translations
      # are not endpoint-preserving diffeomorphisms. The label was wrong, not
      # the arithmetic.
      #
      # What WAS wrong: the boundary rule. A shift estimated by CIRCULAR
      # cross-correlation (the periodic branch above) was then applied with
      # approx(rule = 2), i.e. constant extrapolation -- so a curve shifted by
      # s had its first |s| of the domain filled with a repeat of the endpoint
      # value, flattening exactly the region the circular estimate said should
      # wrap round from the other end. On 24-hour data with a 2.4 h shift that
      # is a tenth of the cycle replaced by a constant. When the design is
      # periodic the warp now wraps; when it is not, it still clamps, and the
      # amount of extrapolation is reported rather than left implicit.
      s_max <- 0.25   # a shift beyond a quarter of the domain is not identified
      shifts[i] <- max(-s_max, min(s_max, shifts[i]))
      h <- time_points - shifts[i]

      if(periodic) {
        # Wrap into the observed domain: the curve is one period of a cycle.
        span <- diff(range(time_points))
        h_use <- min(time_points) +
                 ((h - min(time_points)) %% span)
        extrap_frac[i] <- 0
      } else {
        h_use <- h
        rng <- range(time_points)
        extrap_frac[i] <- mean(h < rng[1] | h > rng[2])
      }
      warp_functions[,i] <- h_use

      # Apply warping
      if(abs(shifts[i]) > 1e-8) {
        # Interpolate curve at warped time points. rule = 2 is a no-op in the
        # periodic branch, where h_use is inside the domain by construction.
        registered_curves[,i] <- approx(time_points, curves[,i],
                                        xout = h_use,
                                        rule = 2)$y
      } else {
        registered_curves[,i] <- curves[,i]
      }
    }
    
    # Create fd objects for output
    basis <- fd_obj$basis
    reg_smooth <- smooth.basis(time_points, registered_curves, basis)
    
    # Create warping function fd objects
    warp_basis <- create.bspline.basis(rangeval = c(0, 1), nbasis = 10)
    warp_smooth <- smooth.basis(time_points, warp_functions, warp_basis)
    
    return(list(
      regfd = reg_smooth$fd,
      registered_curves = registered_curves,
      warp_functions = warp_functions,
      shifts = shifts,
      extrap_frac = extrap_frac,
      warp_direction = "registered -> original",   # P5.6, module-wide
      boundary = if(periodic) "periodic wrap" else "constant extrapolation",
      method = "linear_shift",
      time_points = time_points
    ))
    
  }, error = function(e) {
    cat("Error in linear_shift_alignment:", e$message, "\n")
    # Return identity warping as fallback
    n_curves <- ncol(fd_obj$coefs)
    n_time <- length(time_points)
    
    return(list(
      regfd = fd_obj,
      registered_curves = eval.fd(time_points, fd_obj),
      warp_functions = matrix(rep(time_points, n_curves), n_time, n_curves),
      shifts = rep(0, n_curves),
      method = "identity",
      time_points = time_points
    ))
  })
}

# Parametric alignment function
parametric_alignment <- function(fd_obj, family = "power", 
                                 param_range = c(0.5, 2), 
                                 symmetric = FALSE,
                                 time_points = seq(0, 1, length.out = 100)) {
  
  tryCatch({
    n_curves <- ncol(fd_obj$coefs)
    n_time <- length(time_points)
    
    curves <- eval.fd(time_points, fd_obj)
    mean_curve <- rowMeans(curves)
    
    warp_functions <- matrix(NA, n_time, n_curves)
    registered_curves <- matrix(NA, n_time, n_curves)
    alpha_values <- numeric(n_curves)
    
    # AUDIT (P0.8): the quadratic family alpha*t^2 + (1-alpha)*t has
    # derivative 2*alpha*t + (1-alpha), which is negative near t = 0 whenever
    # alpha > 1. The UI's default search range was c(0.5, 2), so the optimiser
    # could and did return non-monotone maps; clipping them to [0,1] produced
    # a flat leading segment that collapses many new times onto one old time
    # (at alpha = 2, 12 of 24 grid points mapped to time 0). A time warp that
    # is not strictly increasing is not a reparameterisation of time.
    #
    # AUDIT (P4.1): that fix was right about monotonicity and WRONG about
    # where each family's identity lies, and a reviewer caught both halves.
    #
    #   exponential  h(t) = (e^(at) - 1)/(e^a - 1).  Its identity is the LIMIT
    #                a -> 0 (L'Hopital: (at + O(a^2))/(a + O(a^2)) -> t). The
    #                code special-cased  abs(alpha - 1) < 0.001  and returned
    #                t there. At a = 1 the function is (e^t - 1)/(e - 1),
    #                which is not t -- it is the most-curved member of the
    #                family the default range could reach. So the guard put a
    #                DISCONTINUITY in the objective in the middle of the
    #                default search interval [0.5, 2], and optimize() can and
    #                does converge onto it: a curve reported as "alpha = 1,
    #                identity" had in fact been left unwarped by accident
    #                while its neighbours were warped by a real map. And the
    #                genuinely singular point, a = 0 (0/0), had no guard at
    #                all -- it was simply outside the range the UI allowed.
    #   quadratic    identity at a = 0. The P0.8 clamp forced a >= 0.05 and
    #                the UI default gave [0.5, 1], so the identity was NOT
    #                REACHABLE: a curve needing no registration was deformed
    #                anyway, by at least a = 0.5. The 0.05 floor was copied
    #                from the power family, where it is correct; here it is
    #                not. The family is monotone on the whole of (-1, 1).
    #   logistic     identity as its steepness -> 0, also excluded, and also
    #                0/0 at exactly 0 (L1 - L0 = 0). Monotone for every a != 0.
    #   power        h(t) = t^a, identity at a = 1, monotone for a > 0. The
    #                only family whose identity the old range contained.
    #
    # Each family now declares its identity and the open interval on which it
    # is a strictly increasing bijection of [0,1]. The user's range is clamped
    # to that interval AND widened, if necessary, to contain the identity --
    # so "no warping needed" is always inside the search space, in every
    # family. Where the identity is a removable singularity it is returned
    # exactly rather than divided out.
    fam_spec <- switch(family,
      "power"       = list(identity = 1, lo = 0.05, hi = Inf),
      "exponential" = list(identity = 0, lo = -Inf, hi = Inf),
      "quadratic"   = list(identity = 0, lo = -0.999, hi = 0.999),
      "logistic"    = list(identity = 0, lo = -Inf, hi = Inf),
      list(identity = 0, lo = -Inf, hi = Inf))

    param_range <- sort(as.numeric(param_range))
    param_range <- c(max(param_range[1], fam_spec$lo),
                     min(param_range[2], fam_spec$hi))
    # always able to say "this curve needs no warping"
    param_range <- c(min(param_range[1], fam_spec$identity),
                     max(param_range[2], fam_spec$identity))
    param_range <- c(max(param_range[1], fam_spec$lo),
                     min(param_range[2], fam_spec$hi))
    if (param_range[1] >= param_range[2])
      param_range <- c(fam_spec$identity - 1e-3, fam_spec$identity + 1e-3)

    # Define warping function. Every branch returns a map with h(0) = 0,
    # h(1) = 1, strictly increasing on the declared interval, so no clipping
    # to [0,1] is needed -- the pmin/pmax that used to be here only ever
    # masked a non-monotone map instead of rejecting it. They are kept as a
    # floating-point tidy-up (widths of order 1e-16), not as a repair.
    eps_id <- 1e-8
    warp_func <- function(t, alpha) {
      h <- switch(family,
             "power" = t^alpha,
             "exponential" = {
               # identity is the limit a -> 0, where the expression is 0/0
               if (abs(alpha) < eps_id) t
               else (exp(alpha * t) - 1) / (exp(alpha) - 1)
             },
             "quadratic" = alpha * t^2 + (1 - alpha) * t,
             "logistic" = {
               # identity is the limit a -> 0, where L1 - L0 is 0
               if (abs(alpha) < eps_id) t
               else {
                 L  <- function(x) 1 / (1 + exp(-alpha * (x - 0.5)))
                 L0 <- L(0); L1 <- L(1)
                 (L(t) - L0) / (L1 - L0)
               }
             },
             t)
      pmin(1, pmax(0, h))
    }
    
    # AUDIT (P5.11, found by tests/registration_effectiveness_test.R): the
    # search was a bare optimize() over the whole parameter range.
    # optimize() is golden-section plus parabolic interpolation and it
    # assumes the objective is UNIMODAL on the interval. The registration
    # objective is not. On a sharply peaked curve -- which is what a
    # circadian profile is -- the alignment SSE has a deep, narrow well at
    # the correct parameter surrounded by a wide plateau where the curve has
    # been pushed off its own peak. Measured on an already-aligned sample
    # with the power family: SSE is 0.008 at alpha = 1 and about 20 for
    # alpha anywhere in 0.05-0.5 or 1.5-6, and optimize() on [0.05, 6]
    # returned alpha = 6.000, deforming the time axis by 0.58 of the domain
    # on curves that needed no registration at all.
    #
    # This got WORSE at P4.1, not better: widening the ranges to make each
    # family's identity reachable was right, but a wider interval gives
    # optimize() more plateau to get lost on. Registrations run with the
    # narrow pre-P4.1 ranges were partly protected by luck.
    #
    # A coarse grid scan followed by refinement inside the winning bracket
    # is robust to this and costs about 40 extra evaluations per curve --
    # the same pattern fck_auto_lambda() already uses for the GCV search.
    fck_min_1d <- function(f, lo, hi, n_grid = 41) {
      g <- seq(lo, hi, length.out = n_grid)
      v <- vapply(g, function(z) {
        y <- tryCatch(f(z), error = function(e) NA_real_)
        if (is.finite(y)) y else Inf
      }, numeric(1))
      if (all(!is.finite(v))) return(NA_real_)
      k <- which.min(v)
      a <- g[max(1, k - 1)]; b <- g[min(length(g), k + 1)]
      if (a < b) tryCatch(stats::optimize(f, c(a, b), tol = 1e-5)$minimum,
                          error = function(e) g[k])
      else g[k]
    }

    # Optimize warping for each curve
    for(i in 1:n_curves) {
      # Objective function
      objective <- function(alpha) {
        warped_time <- warp_func(time_points, alpha)
        warped_curve <- approx(time_points, curves[,i], xout = warped_time, rule = 2)$y
        sum((warped_curve - mean_curve)^2, na.rm = TRUE)
      }

      a_hat <- fck_min_1d(objective, param_range[1], param_range[2])
      if (!is.finite(a_hat)) a_hat <- fam_spec$identity
      result <- list(minimum = a_hat)
      alpha_values[i] <- result$minimum
      
      # Apply warping
      warp_functions[,i] <- warp_func(time_points, alpha_values[i])
      registered_curves[,i] <- approx(time_points, curves[,i], 
                                      xout = warp_functions[,i], rule = 2)$y
    }
    
    # Create fd objects
    basis <- fd_obj$basis
    reg_smooth <- smooth.basis(time_points, registered_curves, basis)
    
    return(list(
      regfd = reg_smooth$fd,
      registered_curves = registered_curves,
      warp_functions = warp_functions,
      alpha_values = alpha_values,
      family = family,
      warp_direction = "registered -> original",   # P5.6, module-wide
      param_range_used = param_range,
      identity_alpha = fam_spec$identity,
      method = "parametric",
      time_points = time_points
    ))
    
  }, error = function(e) {
    cat("Error in parametric_alignment:", e$message, "\n")
    # Return identity warping
    n_curves <- ncol(fd_obj$coefs)
    n_time <- length(time_points)
    
    return(list(
      regfd = fd_obj,
      registered_curves = eval.fd(time_points, fd_obj),
      warp_functions = matrix(rep(time_points, n_curves), n_time, n_curves),
      alpha_values = rep(1, n_curves),
      method = "identity",
      time_points = time_points
    ))
  })
}

landmark_alignment_simple <- function(fd_obj, landmarks, time_points,
                                    landmark_points = NULL) {
  tryCatch({
    n_curves <- ncol(fd_obj$coefs)
    n_time <- length(time_points)
    
    curves <- eval.fd(time_points, fd_obj)
    rejection_warning <- NULL
    
    # If landmarks are provided, use them for alignment
    if(!is.null(landmark_points) && nrow(landmark_points) > 0) {
      landmark_times <- landmark_points$x
      n_landmarks <- length(landmark_times)
      
      cat("Using", n_landmarks, "landmarks for alignment\n")
      
      # Find corresponding landmark points in each curve
      warp_functions <- matrix(NA, n_time, n_curves)
      registered_curves <- matrix(NA, n_time, n_curves)
      n_rejected <- 0L; rejected_ids <- integer(0)

      for(i in 1:n_curves) {
        # For each curve, find the actual landmarks (peaks/valleys near the specified times)
        curve_landmarks <- numeric(n_landmarks)
        
        for(j in 1:n_landmarks) {
          # Find local extremum near the landmark time
          search_window <- which(abs(time_points - landmark_times[j]) < 0.1)
          if(length(search_window) > 0) {
            local_values <- curves[search_window, i]
            # Find the peak or valley
            if(j %% 2 == 1) {
              # Look for peak for odd landmarks
              curve_landmarks[j] <- time_points[search_window[which.max(local_values)]]
            } else {
              # Look for valley for even landmarks
              curve_landmarks[j] <- time_points[search_window[which.min(local_values)]]
            }
          } else {
            curve_landmarks[j] <- landmark_times[j]
          }
        }
        
        # P5.5: built and VALIDATED in one place. A subject whose detected
        # landmarks cross does not get a folded warp; it gets the identity,
        # and is counted.
        h <- fck_landmark_warp(landmark_times, curve_landmarks, time_points)
        if (is.null(h)) {
          n_rejected <- n_rejected + 1L
          rejected_ids <- c(rejected_ids, i)
          h <- time_points
        }
        warp_functions[, i] <- h

        # P5.6: the module-wide contract, x_R(t) = x(h(t)).
        registered_curves[, i] <- approx(time_points, curves[, i],
                                         xout = h, rule = 2)$y
      }

      if (n_rejected > 0)
        rejection_warning <- sprintf(
          paste("Landmark registration: %d of %d curve(s) had crossed or duplicated",
                "landmarks, which cannot define a monotone time warp. Those curves",
                "were left UNREGISTERED (identity warp) rather than folded:",
                "subject%s %s. Adjust the landmark positions or use fewer."),
          n_rejected, n_curves, if (n_rejected == 1) "" else "s",
          paste(rejected_ids, collapse = ", "))
    } else {
      # No landmarks provided, use automatic detection
      cat("No manual landmarks provided, using automatic landmark detection\n")
      
      # Simple automatic landmark detection: find common peaks
      mean_curve <- rowMeans(curves)
      
      # Find peaks in mean curve
      peaks <- which(diff(sign(diff(mean_curve))) == -2) + 1
      valleys <- which(diff(sign(diff(mean_curve))) == 2) + 1
      
      # Select up to 3 most prominent landmarks
      if(length(peaks) > 0 || length(valleys) > 0) {
        all_extrema <- sort(c(peaks, valleys))
        if(length(all_extrema) > 3) {
          # Select based on prominence
          prominence <- abs(mean_curve[all_extrema] - mean(mean_curve))
          all_extrema <- all_extrema[order(prominence, decreasing = TRUE)[1:3]]
        }
        landmark_times <- time_points[all_extrema]
      } else {
        # Default landmarks at quartiles
        landmark_times <- c(0.25, 0.5, 0.75)
      }
      
      # Apply landmark registration
      warp_functions <- matrix(NA, n_time, n_curves)
      registered_curves <- matrix(NA, n_time, n_curves)
      
      # ---- AUDIT (P0.8): this branch used to be ---------------------------
      #     distortion <- sin(2*pi*time_points) * runif(1, -0.02, 0.02)
      #     warp_functions[,i]    <- pmin(1, pmax(0, time_points + distortion))
      #     registered_curves[,i] <- curves[,i]
      # under the comment "# Simple landmark alignment". It computed landmark
      # positions, discarded them, returned a warp made entirely of random
      # numbers, and did not warp the curves at all -- while the UI announced
      # "Time-warped PCA completed!".
      #
      # A landmark warp is the monotone piecewise-linear map that carries each
      # curve's own landmarks onto the reference landmarks. Detecting the
      # landmarks is done by the same rule for every curve (the largest local
      # extremum nearest each reference position), so the result is
      # deterministic.
      ref_lm <- sort(unique(pmin(1 - 1e-6, pmax(1e-6, landmark_times))))
      find_lm <- function(v, target) {
        # search a window around the reference position for the strongest
        # turning point; fall back to the reference position itself
        w <- which(abs(time_points - target) <= 0.12)
        if (length(w) < 3) return(target)
        d <- diff(v[w])
        turn <- which(d[-length(d)] * d[-1] < 0)
        if (!length(turn)) return(target)
        best <- turn[which.max(abs(d[turn]))]
        time_points[w[best + 1]]
      }
      n_rejected <- 0L; rejected_ids <- integer(0)
      for(i in 1:n_curves) {
        own <- vapply(ref_lm, function(tt) find_lm(curves[, i], tt), numeric(1))
        # P5.6: knots are (reference -> own), i.e. h maps REGISTERED time to
        # ORIGINAL time -- the same direction as the manual branch. This
        # branch used to build (own -> reference), the inverse map, and then
        # apply it as approx(h, curve, xout = t), which inverts it again. The
        # registered curves were right; warp_functions held the wrong object,
        # and every statistic computed from it was not comparable with the
        # other methods'. P5.5: crossed landmarks are rejected, not clipped.
        h <- fck_landmark_warp(ref_lm, own, time_points)
        if (is.null(h)) {
          n_rejected <- n_rejected + 1L
          rejected_ids <- c(rejected_ids, i)
          h <- time_points
        }
        warp_functions[, i] <- h
        registered_curves[, i] <- approx(time_points, curves[, i],
                                         xout = h, rule = 2)$y
      }

      if (n_rejected > 0)
        rejection_warning <- sprintf(
          paste("Automatic landmark registration: %d of %d curve(s) gave crossed or",
                "duplicated landmarks and were left UNREGISTERED (identity warp)",
                "rather than folded: subject%s %s."),
          n_rejected, n_curves, if (n_rejected == 1) "" else "s",
          paste(rejected_ids, collapse = ", "))
    }
    
    basis <- fd_obj$basis
    reg_smooth <- smooth.basis(time_points, registered_curves, basis)
    
    return(list(
      regfd = reg_smooth$fd,
      registered_curves = registered_curves,
      warp_functions = warp_functions,
      method = "landmark",
      n_rejected = n_rejected,
      rejected_ids = rejected_ids,
      rejection_warning = rejection_warning,
      warp_direction = "registered -> original",   # P5.6, module-wide
      time_points = time_points,
      landmarks_used = if(!is.null(landmark_points) && nrow(landmark_points) > 0)
        landmark_points$x else NULL
    ))
    
  }, error = function(e) {
    cat("Error in landmark_alignment:", e$message, "\n")
    return(NULL)
  })
}
