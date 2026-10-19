# ==============================================================================
# server/05b_helpers_warpstats.R -- how much warping, and is it needed?
#
# A registration criterion always improves as the warps get more flexible, so
# there is no GCV for "how much to warp" (Marron, Ramsay, Sangalli & Srivastava
# 2015, Statistical Science). Three things CAN be said, and are computed here:
#
#   1. dance_shift_model(): whether warping is needed at all, by likelihood.
#      A shape-invariant (self-modelling) mixed model -- every participant's
#      curve is one common curve, moved up or down and SHIFTED in time:
#          y_ij = b0 + a_i + f(t_ij - s_i) + e_ij,
#          a_i ~ N(0, sd_a^2),  s_i ~ N(0, sd_s^2),  e_ij ~ N(0, sigma^2)
#      (Lawton, Sylvestre & Maggio 1972; Kneip & Gasser 1988; Brumback &
#      Lindstrom 2004), f a natural cubic spline, fitted by nlme (ML). The
#      spread of the shifts, sd_s, is ESTIMATED rather than tuned, and the
#      likelihood ratio against the same model without shifts (linear, so
#      fitted exactly by lme) tests whether
#      timing differs between participants at all (sd_s = 0 is on the boundary,
#      so p is half the chi-square(1) tail; Self & Liang 1987).
#   2. dance_amp_phase(): how a registration split the variation, by Kneip &
#      Ramsay's (2008, JASA) amplitude/phase decomposition (fda::AmpPhaseDecomp)
#      -- for warps that map the recording onto itself (parametric, landmark,
#      a periodic shift); not for a shift of a non-periodic recording.
#   3. dance_warp_estimates() / dance_warp_group_test(): each participant's warp
#      as a number, in the time units of the data, and whether groups differ.
# ==============================================================================

# ---- 1. is warping needed? ---------------------------------------------------
# Y: participants x time points (NA allowed); t: the times, in the data's units;
# df: flexibility of the common curve f (natural spline degrees of freedom).
dance_shift_model <- function(Y, t, df = 6, ids = rownames(Y)) {
  if (!requireNamespace("nlme", quietly = TRUE))
    return(list(ok = FALSE, message = "The shift model needs the nlme package."))
  Y <- as.matrix(Y); t <- as.numeric(t)
  if (is.null(ids)) ids <- sprintf("S%03d", seq_len(nrow(Y)))
  d <- data.frame(subject = factor(rep(ids, each = ncol(Y)), levels = unique(ids)),
                  t = rep(t, nrow(Y)), y = as.numeric(t(Y)))
  d <- d[is.finite(d$y), , drop = FALSE]
  df <- as.integer(df)
  if (nlevels(droplevels(d$subject)) < 3L || length(unique(d$t)) < df + 2L)
    return(list(ok = FALSE, message = sprintf(paste(
      "The shift model needs at least 3 participants and %d distinct time points",
      "for a common curve with %d degrees of freedom."), df + 2L, df)))
  d$subject <- droplevels(d$subject)
  nsb <- splines::ns(d$t, df = df)
  # f(x) for any x: a natural spline continues as a straight line beyond the
  # data, which is where a shifted time can land at either end
  fsim <- function(x, B) rowSums(stats::predict(nsb, x) * B)
  bn <- paste0("b", seq_len(df))
  st <- stats::coef(stats::lm(d$y ~ nsb)); names(st) <- c("b0", bn)
  fixed <- stats::as.formula(paste("b0 +", paste(bn, collapse = " + "), "~ 1"))
  B <- paste0("cbind(", paste(bn, collapse = ", "), ")")
  # nlme evaluates the model expression in the data, then in its own frame and
  # the search path -- not in this function -- so the curve function is made
  # visible there for the length of the fit, under a name nothing else uses,
  # and removed again on exit
  fname <- ".dance_shift_model_f"
  had <- exists(fname, envir = globalenv(), inherits = FALSE)
  old_f <- if (had) get(fname, envir = globalenv()) else NULL
  assign(fname, fsim, envir = globalenv())
  on.exit(if (had) assign(fname, old_f, envir = globalenv())
          else rm(list = fname, envir = globalenv()), add = TRUE)
  form1 <- stats::as.formula(paste0("y ~ b0 + ", fname, "(t - s, ", B, ")"))
  ctl <- nlme::nlmeControl(maxIter = 200, pnlsMaxIter = 20, msMaxIter = 200)
  fit <- function(form, random) tryCatch(
    nlme::nlme(form, data = d, fixed = fixed, random = random, groups = ~ subject,
               start = st, method = "ML", control = ctl),
    error = function(e) conditionMessage(e))
  # without shifts the model is linear in its parameters: fitted exactly by lme
  d$X <- stats::predict(nsb, d$t)
  m0 <- tryCatch(nlme::lme(y ~ X, random = ~ 1 | subject, data = d, method = "ML",
                           control = nlme::lmeControl(maxIter = 200, msMaxIter = 200)),
                 error = function(e) conditionMessage(e))
  m1 <- fit(form1, nlme::pdDiag(b0 + s ~ 1))
  if (is.character(m1))
    return(list(ok = FALSE, message = paste(
      "The model with time shifts did not fit:", m1, "-- this usually means the spread",
      "of the shifts is at or near zero (no timing variation to estimate); try a",
      "different spline df to check.")))
  if (is.character(m0))
    return(list(ok = FALSE, message = paste("The model without time shifts did not fit:", m0)))
  ll1 <- as.numeric(stats::logLik(m1)); ll0 <- as.numeric(stats::logLik(m0))
  stat <- max(0, 2 * (ll1 - ll0))
  vc <- suppressWarnings(as.numeric(nlme::VarCorr(m1)[, "StdDev"]))
  names(vc) <- rownames(nlme::VarCorr(m1))
  ci <- tryCatch({
    iv <- nlme::intervals(m1, which = "var-cov")$reStruct$subject
    as.numeric(iv["sd(s)", c("lower", "upper")])
  }, error = function(e) c(NA_real_, NA_real_))
  re <- nlme::ranef(m1)
  list(ok = TRUE, df = df, n_participants = nlevels(d$subject), n_obs = nrow(d),
       sd_shift = unname(vc["s"]), sd_shift_ci = ci, sd_level = unname(vc["b0"]),
       sigma = m1$sigma,
       shifts = data.frame(participant = rownames(re), shift = re[, "s"],
                           stringsAsFactors = FALSE),
       lrt = list(statistic = stat, df = 1L,
                  p = if (stat > 1e-8) 0.5 * stats::pchisq(stat, 1, lower.tail = FALSE) else 1),
       logLik = c(with_shift = ll1, without_shift = ll0),
       AIC = c(with_shift = stats::AIC(m1), without_shift = stats::AIC(m0)),
       BIC = c(with_shift = stats::BIC(m1), without_shift = stats::BIC(m0)),
       time_range = range(d$t), model = m1)
}

dance_shift_model_text <- function(sm, unit = "h") {
  if (!isTRUE(sm$ok)) return(sm$message %||% "The shift model was not fitted.")
  f <- function(x, k = 2) formatC(x, digits = k, format = "f")
  paste0(
    "Shape-invariant model: y_ij = b0 + a_i + f(t_ij - s_i) + e_ij\n",
    "  one common curve f (natural cubic spline, ", sm$df, " df), moved up or down (a_i)\n",
    "  and shifted in time (s_i) for each participant; fitted by maximum likelihood\n",
    "  (nlme) to ", sm$n_obs, " observations of ", sm$n_participants, " participants.\n\n",
    "Spread of the time shifts between participants:\n",
    "  SD(s) = ", f(sm$sd_shift), " ", unit,
    if (all(is.finite(sm$sd_shift_ci)))
      paste0("  (95% CI ", f(sm$sd_shift_ci[1]), " to ", f(sm$sd_shift_ci[2]), " ", unit, ")") else "",
    "\n\nIs warping needed? Shifts against none (likelihood ratio):\n",
    "  chi2(1) = ", f(sm$lrt$statistic), ", p ",
    if (sm$lrt$p < 0.001) "< .001" else paste("=", formatC(sm$lrt$p, digits = 3, format = "f")),
    "  (half the chi-square(1) tail: SD(s) = 0 is on the boundary)\n",
    "  AIC ", f(sm$AIC["with_shift"], 1), " with shifts vs ", f(sm$AIC["without_shift"], 1), " without;",
    " BIC ", f(sm$BIC["with_shift"], 1), " vs ", f(sm$BIC["without_shift"], 1), "\n\n",
    "Reading it: a significant ratio says participants' curves differ in TIMING,\n",
    "beyond their level, so a shift registration has something real to align. SD(s)\n",
    "says how much: a registration moving curves far beyond about 2 x SD(s) is\n",
    "removing more than timing. A non-significant ratio says the timing differences\n",
    "are not distinguishable from noise around one curve, and warping is not\n",
    "supported. The model allows a shift only (no stretching of time) and treats\n",
    "residuals as independent. The predicted shifts per participant are shrunk\n",
    "towards zero (they are random-effect predictions).")
}

# ---- 2. how a registration split the variation -------------------------------
# xfd: the curves before registration; warp: a registration result (regfd,
# warp_functions, time_points, method). For a translation h(t) = t - s the
# derivative is 1 everywhere, so the identity is passed as h (fda reads only its
# derivative); a wrapped periodic h would otherwise show a jump at the wrap.
dance_amp_phase <- function(xfd, warp) {
  if (is.null(warp$regfd) || is.null(warp$time_points))
    return(list(ok = FALSE, message = "No registered curves."))
  tp <- warp$time_points
  translation <- warp$method %in% c("linear_shift", "identity")
  # The decomposition rests on a change of variables that holds only for warps
  # mapping the recording onto itself. A shift on a cycle does (it wraps); a
  # shift of a non-periodic recording slides past its ends, and the result is
  # not a decomposition (on a test it put -107% of the variation in phase).
  if (translation && !identical(warp$boundary, "periodic wrap") && !identical(warp$method, "identity"))
    return(list(ok = FALSE, message = paste(
      "not defined for a shift of a non-periodic recording: a translation slides the",
      "curves past the ends of the recording instead of mapping it onto itself, which",
      "the decomposition requires. Use the shift model's SD(s) and the per-participant",
      "shifts instead.")))
  H <- if (translation) matrix(tp, length(tp), ncol(xfd$coefs))
       else as.matrix(warp$warp_functions)
  hfd <- dance_refd(H, tp, NULL)
  rng <- range(tp)
  out <- tryCatch(fda::AmpPhaseDecomp(xfd, warp$regfd, hfd, rng),
                  error = function(e) conditionMessage(e))
  if (is.character(out)) return(list(ok = FALSE, message = out))
  c(list(ok = TRUE, translation = translation), out)
}

dance_amp_phase_text <- function(ap) {
  if (!isTRUE(ap$ok)) return(paste("Amplitude/phase decomposition not available:", ap$message %||% ""))
  paste0(
    "Amplitude/phase decomposition (Kneip & Ramsay 2008; fda::AmpPhaseDecomp):\n",
    "  MS amplitude = ", formatC(ap$MS.amp, digits = 4, format = "g"),
    "   MS phase = ", formatC(ap$MS.pha, digits = 4, format = "g"),
    "   C = ", formatC(ap$C, digits = 4, format = "f"), "\n",
    "  Share of the variation this registration attributes to phase: ",
    sprintf("%.1f%%", 100 * ap$RSQR), "\n",
    if (isTRUE(ap$translation)) "  (a translation: its derivative is 1, so C is exactly 1)\n" else "",
    "  A large phase share from warps that look implausible (see the shift model and\n",
    "  the per-participant warps) means the registration is absorbing amplitude.\n",
    if (is.finite(ap$RSQR) && ap$RSQR < 0)
      "  A NEGATIVE share: registration increased the mean-square variation.\n" else "")
}

# ---- 3. each participant's warp as a number, and group differences -----------
# span: the length of the time axis the [0, 1] warping domain stands for, in
# `unit`; ids: the curves' participant ids.
dance_warp_estimates <- function(warp, ids, span = 1, unit = "of the recording") {
  tp <- warp$time_points
  H <- as.matrix(warp$warp_functions)
  # signed mean displacement of the time axis: positive = the participant's
  # curve was read LATER in time to align (it ran late)
  disp <- if (!is.null(H) && length(tp) > 1) {
    dt <- diff(tp)
    apply(H - tp, 2, function(v) sum((v[-1] + v[-length(v)]) / 2 * dt) / diff(range(tp)))
  } else rep(NA_real_, length(ids))
  m <- warp$method %||% "identity"
  if (m %in% c("linear_shift", "identity") && !is.null(warp$shifts)) {
    # h(t) = t - s: the curve is read at t - s, so it ran -s late
    est <- -warp$shifts * span
    out <- data.frame(participant = ids, estimate = est, stringsAsFactors = FALSE)
    attr(out, "label") <- sprintf("time shift (%s; positive = later than the reference)", unit)
    attr(out, "circular") <- identical(warp$boundary, "periodic wrap")
  } else if (identical(m, "parametric") && !is.null(warp$alpha_values)) {
    out <- data.frame(participant = ids, estimate = warp$alpha_values,
                      mean_displacement = disp * span, stringsAsFactors = FALSE)
    attr(out, "label") <- sprintf("warp parameter alpha (%s family; identity at %s)",
                                  warp$family %||% "", format(warp$identity_alpha %||% NA))
    attr(out, "circular") <- FALSE
  } else {
    out <- data.frame(participant = ids, estimate = disp * span, stringsAsFactors = FALSE)
    attr(out, "label") <- sprintf("mean displacement of the time axis (%s; positive = later)", unit)
    attr(out, "circular") <- FALSE
  }
  attr(out, "span") <- span
  out
}

# Between-group comparison of the per-participant warps. Linear: Welch's one-way
# ANOVA (stats::oneway.test) and Kruskal-Wallis (stats::kruskal.test). A
# periodic shift is an angle on the cycle: Watson-Williams test.
dance_warp_group_test <- function(est, groups, circular = FALSE, cycle = 1) {
  ok <- is.finite(est) & !is.na(groups)
  est <- est[ok]; g <- droplevels(factor(groups[ok]))
  if (nlevels(g) < 2L) return(list(ok = FALSE, message = "Fewer than two groups."))
  if (any(table(g) < 2L)) return(list(ok = FALSE, message = "Every group needs at least two participants."))
  desc <- do.call(rbind, lapply(levels(g), function(l) {
    x <- est[g == l]
    data.frame(group = l, n = length(x), mean = mean(x), sd = stats::sd(x), median = stats::median(x),
               stringsAsFactors = FALSE)
  }))
  if (isTRUE(circular)) {
    ang <- 2 * pi * est / cycle
    ww <- dance_watson_williams_test(split(ang, g))
    return(list(ok = TRUE, circular = TRUE, groups = desc,
                tests = data.frame(test = "Watson-Williams (circular)", statistic = ww$F,
                                   df = sprintf("%s, %s", ww$df1, ww$df2), p = ww$p,
                                   stringsAsFactors = FALSE)))
  }
  w <- stats::oneway.test(est ~ g, var.equal = FALSE)
  k <- stats::kruskal.test(est ~ g)
  list(ok = TRUE, circular = FALSE, groups = desc,
       tests = data.frame(
         test = c("Welch one-way ANOVA", "Kruskal-Wallis"),
         statistic = c(unname(w$statistic), unname(k$statistic)),
         df = c(sprintf("%.0f, %.1f", w$parameter[1], w$parameter[2]), sprintf("%d", as.integer(k$parameter))),
         p = c(w$p.value, k$p.value), stringsAsFactors = FALSE))
}
