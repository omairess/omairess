# FaRiMa 1.2 — corrected mathematical core, original preliminary interface retained.
# Based on Farimato(prelim).R supplied on 18 September 2026.
# Published model defaults are distinguished from unvalidated composite choices.
# options(farimato.launch=FALSE) loads definitions without launching the app.
required_pkgs <- c("shiny","shinydashboard","shinyWidgets","DT","plotly","ggplot2","dplyr","tidyr","readr","readxl","haven","lubridate","RColorBrewer","deSolve","digest")
missing <- required_pkgs[!vapply(required_pkgs,requireNamespace,logical(1),quietly=TRUE)]
if(length(missing))stop("Install required packages: ",paste(missing,collapse=", "),call.=FALSE)
suppressPackageStartupMessages({library(shiny);library(shinydashboard);library(shinyWidgets);library(DT);library(plotly);library(dplyr)})

# Verified published-model core — September 2026
# Published McCauley model variants; author-defined composite impairment index.
# Run this file with Rscript, or use shiny::runApp() on the supplied app folder.
# To load functions without launching: options(farimato.launch = FALSE); source(...)
# Sources and assumptions are documented beside the relevant controls and in exports.

FR_VERSION <- "1.2-original-interface"
FR_REFERENCES <- c(
  model_2021 = "https://doi.org/10.1016/j.trf.2021.04.006",
  model_2024jtb = "https://doi.org/10.1016/j.jtbi.2024.111851",
  model_2024fe = "https://doi.org/10.3389/fenvh.2024.1362755",
  inherited_PVT = "https://doi.org/10.5665/sleep.3246",
  lineage = "https://doi.org/10.1016/j.jtbi.2008.09.012",
  threshold_exposure = "https://doi.org/10.3357/ASEM.3455.2013",
  KSS = "https://doi.org/10.3109/00207459008994241",
  PVT = "https://doi.org/10.3758/BF03200977",
  solver = "https://doi.org/10.18637/jss.v033.i09")

# One-sentence, plain-language meaning of each published parameter (per the
# model papers in FR_REFERENCES), shown next to the sliders/values so users
# don't have to cross-reference the paper to know what they're changing.
FR_PARAM_INFO <- c(
  p_f = "Floor of the outcome scale: the score predicted with no sleep pressure, no sleep inertia and the circadian drive at its least alerting point.",
  alpha_w = "Rate constant (1/h) at which the homeostatic process climbs toward its wake asymptote the longer someone stays awake.",
  alpha_s = "Rate constant (1/h) at which the homeostatic process falls toward its sleep asymptote the longer someone stays asleep.",
  beta = "Coupling strength from the allostatic (multi-day sleep debt) process into the homeostatic process's equilibrium point.",
  W_c = "Critical wake duration (h) used to derive the sleep recovery rate of the allostatic process from its wake accumulation rate.",
  eta_w = "Rate constant (1/h) at which the allostatic process (slow, multi-day sleep debt) builds up during wake.",
  zeta = "Strength of the reciprocal feedback from the homeostatic process back onto the allostatic process; only present in the 2024 (Frontiers) model version, 0 otherwise.",
  phi = "Clock time (h) of the circadian sine's peak phase, i.e. roughly when the circadian drive is most alerting.",
  lambda = "Rate constant (1/h) at which the circadian amplitude relaxes toward full strength after a sleep/wake transition.",
  mu_s = "Circadian offset added to the sine drive while asleep; shifts how the circadian rhythm modulates the homeostatic process overnight.",
  mu_w = "Circadian offset added to the sine drive while awake; shifts how the circadian rhythm modulates the homeostatic process during the day.",
  nu = "Decay rate (1/h) of sleep inertia after waking, i.e. how fast post-sleep grogginess wears off.",
  gamma = "How strongly sleep inertia feeds back into (slows) homeostatic recovery while asleep.",
  xi = "Scaling factor that converts the circadian/inertia drive into the outcome metric's own units; calibrated separately for KSS and PVT.")

# Parameters offered to the "Suggest calibration" optimizer, before per-tab
# filtering by model version (gamma/nu absent for "2021"; zeta absent unless
# "2024fe" -- see make_calibration_suggester).
FR_CALIB_PARAM_CHOICES <- c("alpha_w","alpha_s","W_c","mu_w","p_f","xi",
                             "eta_w","beta","zeta","phi","lambda","mu_s","nu","gamma")

fr_scalar <- function(x, name, lo = -Inf, hi = Inf, open_lo = FALSE, open_hi = FALSE) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) ||
      if (open_lo) x <= lo else x < lo) stop(name, " is outside its valid range.", call. = FALSE)
  if (if (open_hi) x >= hi else x > hi) stop(name, " is outside its valid range.", call. = FALSE)
  x
}

`%||%` <- function(a, b) if (is.null(a)) b else a

mccauley_params <- function(metric = c("PVT", "KSS"), version = c("2024fe", "2024jtb", "2021")) {
  metric <- match.arg(metric); version <- match.arg(version); pvt <- metric == "PVT"
  P <- list(metric = metric, version = version, p_f = if (pvt) 0 else 1,
    alpha_w = if (pvt) .028 else .22, alpha_s = if (pvt) .26 else .037,
    W_c = if (pvt) 20.2 else 22.02, xi = if (pvt) 1.09 else .51,
    beta = .26, mu_s = -1.5, phi = 21.2, lambda = .49, tau = 24, T_day = 24,
    eta_w = if (version == "2024fe") .0126 else .0074,
    mu_w = if (version == "2024fe") { if (pvt) .466 else .82 } else .33,
    zeta = if (version == "2024fe") 1.31 else 0,
    nu = 1.37, gamma = if (version == "2021") 0 else .71,
    use_inertia = version != "2021")
  P$eta_s <- P$eta_w * P$W_c / (24 - P$W_c)
  P
}

# Proportional exploration band around a study's published value, used to seed
# the "Explore custom parameters" sliders. Bounded parameters (phi, W_c) get an
# additive band clamped to their physical domain instead of a multiplicative one.
fr_slider_range <- function(nm, value) {
  bounded <- list(phi = c(0, 24), W_c = c(1e-3, 24 - 1e-3))
  if (nm %in% names(bounded)) {
    b <- bounded[[nm]]; span <- (b[2] - b[1]) * 0.35
    return(c(max(b[1], value - span), min(b[2], value + span)))
  }
  if (nm %in% c("beta", "gamma", "zeta")) return(c(0, if (value == 0) 1 else value * 3))
  if (nm %in% c("alpha_w", "alpha_s", "eta_w", "lambda", "xi", "nu")) return(c(value * 0.2, value * 3))
  half <- if (value == 0) 1 else abs(value)
  c(value - half, value + half)
}

fr_parameters <- function(metric, version, overrides = list()) {
  P <- mccauley_params(metric, version)
  allowed <- c("p_f", "alpha_w", "alpha_s", "W_c", "xi", "beta", "mu_s", "phi",
               "lambda", "eta_w", "mu_w", "zeta", "nu", "gamma")
  if (length(setdiff(names(overrides), allowed))) stop("Unknown parameter override.", call. = FALSE)
  for (nm in names(overrides)) P[[nm]] <- overrides[[nm]]
  # Structural exclusions belong to the selected model, even after custom edits.
  if (version != "2024fe") P$zeta <- 0
  if (version == "2021") { P$gamma <- 0; P$nu <- 1.37; P$use_inertia <- FALSE }
  for (nm in c("alpha_w", "alpha_s", "eta_w", "lambda", "xi", "nu"))
    fr_scalar(P[[nm]], nm, 0, open_lo = TRUE)
  for (nm in c("beta", "gamma", "zeta")) fr_scalar(P[[nm]], nm, 0)
  for (nm in c("p_f", "mu_s", "mu_w")) fr_scalar(P[[nm]], nm)
  fr_scalar(P$phi, "phi", 0, 24, open_hi = TRUE)
  fr_scalar(P$W_c, "W_c", 0, 24, open_lo = TRUE, open_hi = TRUE)
  P$eta_s <- P$eta_w * P$W_c / (24 - P$W_c)
  fr_scalar(P$eta_s, "eta_s", 0, open_lo = TRUE)
  P
}

fr_empty_intervals <- function() data.frame(onset_h = numeric(), offset_h = numeric())

fr_intervals <- function(x, end, name = "Sleep schedule") {
  fr_scalar(end, "Study duration", 0, open_lo = TRUE)
  if (!is.data.frame(x) || !all(c("onset_h", "offset_h") %in% names(x)))
    stop(name, " needs onset_h and offset_h columns.", call. = FALSE)
  x <- x[, c("onset_h", "offset_h"), drop = FALSE]
  if (!all(vapply(x, is.numeric, logical(1))) || any(!is.finite(as.matrix(x))))
    stop(name, " contains nonnumeric or nonfinite boundaries.", call. = FALSE)
  if (any(x$onset_h < 0 | x$offset_h > end | x$offset_h <= x$onset_h))
    stop(name, ": each interval must satisfy 0 <= onset < offset <= study duration.", call. = FALSE)
  x <- x[order(x$onset_h, x$offset_h), , drop = FALSE]; rownames(x) <- NULL
  if (nrow(x) > 1L && any(x$onset_h[-1L] < head(x$offset_h, -1L)))
    stop(name, " has overlapping intervals. Resolve them explicitly.", call. = FALSE)
  x
}

fr_inside <- function(t, intervals) {
  ans <- rep(FALSE, length(t))
  for (i in seq_len(nrow(intervals)))
    ans <- ans | (t >= intervals$onset_h[i] & t < intervals$offset_h[i])
  ans
}

fr_rhs <- function(t, y, P, wake) {
  q <- y[1]; u <- y[2]; h <- y[3]; k <- y[4]
  c <- sin(2*pi*(t - P$phi)/24)
  if (wake) {
    dy <- c(-P$alpha_w*(q-P$beta*u) + P$xi*k*(c+P$mu_w),
            P$eta_w*(u-P$zeta*P$xi*q), -P$nu*h, P$lambda*k*(1-k))
  } else {
    dy <- c(-P$alpha_s*(q-P$beta*(u-1/P$eta_s)) + P$xi*k*(c+P$mu_s),
            P$eta_s*(P$zeta*P$xi*q-u)+1, -P$nu*(h-P$gamma*q), -P$lambda*k)
  }
  if (!P$use_inertia) dy[3] <- 0
  list(unname(dy))
}

fr_advance <- function(y, times, P, wake) {
  stopifnot(length(y) == 4L, all(is.finite(y)), length(times) >= 2L, all(diff(times) > 0))
  y <- setNames(as.numeric(y), c("q", "u", "h", "k"))
  z <- deSolve::ode(y, times, fr_rhs, parms = P, wake = wake, method = "lsoda",
                    rtol = 1e-10, atol = 1e-11, maxsteps = 100000)
  if (nrow(z) != length(times) || any(!is.finite(z)) || attr(z, "istate")[1] < 0)
    stop("The solver did not return a complete finite trajectory.", call. = FALSE)
  if (any(z[, "k"] < -1e-8 | z[, "k"] > 1+1e-8))
    stop("Circadian amplitude left its valid range; check the model settings.", call. = FALSE)
  # Only remove solver roundoff, never impose a positive physiological floor.
  z[, "k"] <- pmin(pmax(z[, "k"], 0), 1)
  z
}

fr_endpoint <- function(y, a, b, P, wake) {
  if (b <= a) return(as.numeric(y))
  as.numeric(tail(fr_advance(y, c(a, b), P, wake), 1)[, -1])
}

fr_baseline <- function(P, onset = 23, offset = 7, t0_clock = 8) {
  fr_scalar(onset, "Baseline onset", 0, 24, open_hi = TRUE)
  fr_scalar(offset, "Baseline awakening", 0, 24, open_hi = TRUE)
  fr_scalar(t0_clock, "Study start clock", 0, 24, open_hi = TRUE)
  S <- (offset-onset) %% 24; W <- 24-S
  if (S == 0) stop("Baseline onset and awakening must differ.", call. = FALSE)
  # Analytic periodic normalized kappa. At W <= S its limiting value is zero.
  k0 <- if (W > S) expm1(P$lambda*(W-S))/expm1(P$lambda*W) else 0
  if (!is.finite(k0)) stop("Circadian rate is too large for baseline calculation.", call. = FALSE)
  day_map <- function(x) {
    y <- c(x, if (length(x) == 2L) 0, k0)
    y <- fr_endpoint(y, offset, offset+W, P, TRUE)
    fr_endpoint(y, offset+W, offset+24, P, FALSE)
  }
  n <- if (P$use_inertia) 3L else 2L
  zero <- rep(0, n); b <- day_map(zero)[seq_len(n)]
  M <- vapply(seq_len(n), function(j) {
    x <- zero; x[j] <- 1
    day_map(x)[seq_len(n)]-b
  }, numeric(n))
  rho <- max(Mod(eigen(M, only.values = TRUE)$values))
  if (!is.finite(rho) || rho >= 1-1e-10)
    stop(sprintf("%s baseline is not a stable periodic state (cycle multiplier %.8f). Use more baseline sleep or different parameters.",
                 P$metric, rho), call. = FALSE)
  x <- tryCatch(as.numeric(solve(diag(n)-M, b)), error = function(e)
    stop("Baseline fixed-point system is singular.", call. = FALSE))
  state <- c(x, if (n == 2L) 0, k0)
  check <- day_map(x)
  residual <- max(abs(check-state))
  tol <- 1e-8 + 1e-9*max(abs(state))
  if (!all(is.finite(state)) || !is.finite(residual) || residual > tol)
    stop("Baseline solution failed its one-day residual check.", call. = FALSE)
  elapsed <- (t0_clock-offset) %% 24
  shifted <- fr_endpoint(state, offset, offset+min(elapsed, W), P, TRUE)
  if (elapsed > W) shifted <- fr_endpoint(shifted, offset+W, offset+elapsed, P, FALSE)
  list(state = shifted, awakening_state = state, residual = residual, tolerance = tol,
       cycle_multiplier = rho, onset = onset, offset = offset,
       method = "Stable affine periodic solution, checked over a complete baseline day")
}

fr_initial_state <- function(P, initial) {
  if (is.null(initial) || !all(c("p", "u", "h", "kappa") %in% names(initial)))
    stop("Supply explicit p, u, h and kappa initial states for both outcomes.", call. = FALSE)
  z <- vapply(c("p", "u", "h", "kappa"), function(nm) fr_scalar(initial[[nm]], nm), numeric(1))
  if (z["kappa"] < 0 || z["kappa"] > 1) stop("Initial kappa must be in [0,1].", call. = FALSE)
  if (!P$use_inertia && z["h"] != 0) stop("The 2021 initial h must be zero.", call. = FALSE)
  unname(c(z["p"]-P$p_f-P$xi*z["h"], z["u"], z["h"], z["kappa"]))
}

fr_simulate_arm <- function(P, schedule, end, t0_clock, spacing = 1/60,
                            baseline = list(onset = 23, offset = 7), initial = NULL,
                            extra_times = numeric()) {
  fr_scalar(spacing, "Output spacing", 1/3600, .25)
  fr_scalar(t0_clock, "Study start clock", 0, 24, open_hi = TRUE)
  schedule <- fr_intervals(schedule, end)
  if (end/spacing > 500000) stop("Too many output rows; shorten the study or increase output spacing.", call. = FALSE)
  if (is.null(baseline)) {
    init <- list(state = fr_initial_state(P, initial), method = "User-supplied initial state")
  } else init <- fr_baseline(P, baseline$onset, baseline$offset, t0_clock)
  bounds <- sort(unique(c(0, end, unlist(schedule, use.names = FALSE))))
  extra_times <- extra_times[is.finite(extra_times) & extra_times >= 0 & extra_times <= end]
  grid <- sort(unique(c(seq(0, end, by = spacing), seq(0, floor(end), by = 1), bounds, extra_times)))
  state <- init$state; pieces <- vector("list", length(bounds)-1L)
  for (i in seq_len(length(bounds)-1L)) {
    a <- bounds[i]; b <- bounds[i+1L]
    wake <- !fr_inside((a+b)/2, schedule)
    tt <- sort(unique(c(a, grid[grid > a & grid < b], b)))
    z <- fr_advance(state, t0_clock+tt, P, wake)
    state <- as.numeric(tail(z, 1)[, -1])
    z[, 1] <- tt
    if (i > 1L) z <- z[-1L, , drop = FALSE]
    pieces[[i]] <- z
  }
  z <- do.call(rbind, pieces)
  d <- data.frame(time_h = z[,1], clock_h = (t0_clock+z[,1]) %% 24,
    wake = !fr_inside(z[,1], schedule), p = P$p_f+z[,2]+P$xi*z[,4],
    q = z[,2], u = z[,3], h = z[,4], kappa = z[,5],
    c = sin(2*pi*(t0_clock+z[,1]-P$phi)/24), row.names = NULL)
  if (any(!is.finite(as.matrix(d)))) stop("Nonfinite model output.", call. = FALSE)
  list(predictions = d, initialization = init)
}

fr_index_settings <- function(settings) {
  w <- settings$weights
  if (!is.numeric(w) || length(w) != 2L || any(!is.finite(w)) || any(w < 0) || sum(w) <= 0)
    stop("Choose two nonnegative index weights with a positive sum.", call. = FALSE)
  fr_scalar(settings$pvt_ceiling, "Index PVT ceiling", 0, open_lo = TRUE)
  scaled <- w/max(w)
  settings$normalized_weights <- scaled/sum(scaled)
  settings
}

fr_index <- function(pvt, kss, settings) {
  w <- settings$normalized_weights
  w[1]*(pmin(pmax(kss, 1), 9)-1)/8 + w[2]*pmin(pmax(pvt, 0)/settings$pvt_ceiling, 1)
}

fr_positive_segment <- function(y0, y1, width, threshold) {
  a <- y0-threshold; b <- y1-threshold
  if (a <= 0 && b <= 0) return(c(hours = 0, area = 0))
  if (a >= 0 && b >= 0) return(c(hours = width, area = width*(a+b)/2))
  fraction <- if (a > 0) a/(a-b) else b/(b-a)
  c(hours = width*fraction, area = width*fraction*max(a,b)/2)
}

fr_exposure <- function(pred, schedule, duty, settings, index_settings) {
  metric <- match.arg(settings$metric, c("Index", "PVT_raw", "KSS_raw"))
  threshold <- fr_scalar(settings$threshold, "Exposure threshold")
  if (metric == "Index") fr_scalar(threshold, "Index threshold", 0, 1)
  if (metric == "PVT_raw") fr_scalar(threshold, "PVT threshold", 0)
  if (metric == "KSS_raw") fr_scalar(threshold, "KSS threshold", 1, 9)
  scope <- match.arg(settings$scope, c("wake", "duty"))
  n <- nrow(pred)-1L; left <- seq_len(n); width <- diff(pred$time_h)
  mid <- (head(pred$time_h,-1)+tail(pred$time_h,-1))/2
  active <- !fr_inside(mid, schedule)
  if (scope == "duty") active <- active & fr_inside(mid, duty)
  ans <- data.frame(start_h = head(pred$time_h,-1), end_h = tail(pred$time_h,-1),
                    included = active, hours = ifelse(active, width, 0),
                    hours_above = 0, excess_area = 0, index_area = 0, kss_raw_area = 0,
                    kss_bounded_area = 0, pvt_raw_area = 0)
  for (i in which(active)) {
    p <- pred$PVT_raw[c(i,i+1L)]; k <- pred$KSS_raw[c(i,i+1L)]
    # Split linear interpolation wherever a component clips, so composite
    # integration remains piecewise linear instead of bridging its corners.
    cuts <- c(0,1)
    if (diff(p) != 0) cuts <- c(cuts, (c(0,index_settings$pvt_ceiling)-p[1])/diff(p))
    if (diff(k) != 0) cuts <- c(cuts, (c(1,9)-k[1])/diff(k))
    cuts <- sort(unique(cuts[cuts >= 0 & cuts <= 1]))
    pp <- p[1]+cuts*diff(p); kk <- k[1]+cuts*diff(k)
    ii <- fr_index(pp, kk, index_settings)
    yy <- switch(metric, Index = ii, PVT_raw = pp, KSS_raw = kk)
    for (j in seq_len(length(cuts)-1L)) {
      h <- width[i]*(cuts[j+1L]-cuts[j])
      z <- fr_positive_segment(yy[j], yy[j+1L], h, threshold)
      ans$hours_above[i] <- ans$hours_above[i]+z[1]
      ans$excess_area[i] <- ans$excess_area[i]+z[2]
      ans$index_area[i] <- ans$index_area[i]+h*(ii[j]+ii[j+1L])/2
      ans$kss_raw_area[i] <- ans$kss_raw_area[i]+h*(kk[j]+kk[j+1L])/2
      ans$pvt_raw_area[i] <- ans$pvt_raw_area[i]+h*(pp[j]+pp[j+1L])/2
      clipped <- pmin(pmax(kk[c(j,j+1L)],1),9)
      ans$kss_bounded_area[i] <- ans$kss_bounded_area[i]+h*mean(clipped)
    }
  }
  H <- sum(ans$hours); A <- sum(ans$excess_area)
  safe_mean <- function(v) if (H > 0) sum(v)/H else NA_real_
  list(intervals = ans, summary = data.frame(metric = metric, threshold = threshold, scope = scope,
    exposure_hours = H, hours_above_threshold = sum(ans$hours_above), excess_area = A,
    mean_excess = if (H > 0) A/H else NA_real_, mean_index = safe_mean(ans$index_area),
    mean_KSS_raw = safe_mean(ans$kss_raw_area), mean_KSS_bounded = safe_mean(ans$kss_bounded_area),
    mean_PVT_raw = safe_mean(ans$pvt_raw_area),
    sampled_peak_PVT_raw = if (H > 0) max(pred$PVT_raw[unique(c(which(active),which(active)+1L))]) else NA_real_,
    area_units = switch(metric, Index = "index-hours", PVT_raw = "lapses/10-min x hours", KSS_raw = "KSS-points x hours")))
}

fr_parse_numbers <- function(x, decimal = ".", name = "Values", allow_missing = TRUE) {
  decimal <- match.arg(decimal, c(".", ","))
  if (is.numeric(x)) out <- as.numeric(x) else {
    s <- trimws(as.character(x)); missing <- is.na(s) | s %in% c("", "NA")
    # Strict grammar: no thousands separators or partial numeric strings.
    mark <- if (decimal == ".") "\\." else ","
    pattern <- paste0("^[+-]?(?:[0-9]+(?:",mark,"[0-9]*)?|",mark,"[0-9]+)(?:[eE][+-]?[0-9]+)?$")
    bad <- !missing & !grepl(pattern, s, perl = TRUE)
    if (any(bad)) stop(name, ": invalid number(s) at row(s) ", paste(head(which(bad),8),collapse=", "),
                       ". Check the decimal mark; grouping separators are not accepted.", call. = FALSE)
    if (decimal == ",") s <- gsub(",", ".", s, fixed = TRUE)
    out <- suppressWarnings(as.numeric(s)); out[missing] <- NA_real_
  }
  if (any(!is.na(out) & !is.finite(out)) || any(is.nan(out))) stop(name, " contains nonfinite values.", call. = FALSE)
  if (!allow_missing && anyNA(out)) stop(name, " contains missing values.", call. = FALSE)
  out
}

fr_parse_times <- function(text) {
  if (!nzchar(trimws(text))) return(numeric())
  fr_parse_numbers(strsplit(trimws(text), "[;[:space:]]+")[[1]], name = "Elapsed hours", allow_missing = FALSE)
}

fr_parse_iso <- function(x, name = "Datetime") {
  s <- trimws(as.character(x))
  pat <- "^[0-9]{4}-[0-9]{2}-[0-9]{2}[T ][0-9]{2}:[0-9]{2}(:[0-9]{2}(\\.[0-9]+)?)?(Z|[+-][0-9]{2}:?[0-9]{2})$"
  if (anyNA(s) || any(!grepl(pat,s,perl=TRUE)))
    stop(name, ": use ISO datetimes with an explicit UTC offset, e.g. 2026-09-16T08:00:00+02:00. Unqualified or ambiguous local times are not accepted.", call. = FALSE)
  hh <- as.numeric(substr(s,12,13)); mm <- as.numeric(substr(s,15,16))
  local <- sub("(Z|[+-][0-9]{2}:?[0-9]{2})$", "", s)
  seconds <- ifelse(nchar(local)>16, suppressWarnings(as.numeric(substring(local,18))), 0)
  if (any(hh>23 | mm>59 | seconds>=60)) stop(name, " has an invalid clock time.",call.=FALSE)
  numeric_offsets <- s[!grepl("Z$",s)]
  if (length(numeric_offsets)) {
    offset <- gsub(":", "", sub("^.*([+-][0-9]{2}:?[0-9]{2})$", "\\1", numeric_offsets))
    oh <- as.numeric(substr(offset,2,3)); om <- as.numeric(substr(offset,4,5))
    if (any(oh>14 | om>59 | (oh==14 & om!=0))) stop(name, " has an invalid UTC offset.",call.=FALSE)
  }
  short <- nchar(local)==16
  s[short] <- paste0(substr(s[short],1,16),":00",substring(s[short],17))
  z <- suppressWarnings(lubridate::ymd_hms(s, tz = "UTC", quiet = TRUE))
  if (anyNA(z)) stop(name, " contains an invalid date/time.", call. = FALSE)
  z
}

fr_origin <- function(mode = "elapsed", clock = 8, start_iso = NULL) {
  if (mode == "dated") {
    instant <- fr_parse_iso(start_iso, "Study start")
    if (length(instant) != 1L) stop("Supply exactly one study start.", call. = FALSE)
    s <- sub("^[^T ]+[T ]", "", start_iso)
    hh <- as.numeric(substr(s,1,2)); mm <- as.numeric(substr(s,4,5))
    sec <- if (substr(s,6,6) == ":") as.numeric(sub("^([0-9]+(?:\\.[0-9]+)?).*$","\\1",substring(s,7),perl=TRUE)) else 0
    clock <- hh+mm/60+sec/3600
    fr_scalar(clock, "Start clock", 0, 24, open_hi = TRUE)
    list(mode = mode, t0_clock = clock, start_iso = start_iso, start_utc = instant,
         clock_policy = "Elapsed UTC hours; biological clock continues from the starting local clock without DST jumps")
  } else {
    fr_scalar(clock, "Start clock", 0, 24, open_hi = TRUE)
    list(mode = "elapsed", t0_clock = clock, start_iso = NULL, start_utc = NULL,
         clock_policy = "Fixed 24-hour model clock; elapsed hours from the explicit study origin")
  }
}

fr_elapsed <- function(x, format, origin, decimal = ".", name = "Time") {
  if (format == "hours") return(fr_parse_numbers(x, decimal, name, FALSE))
  if (is.null(origin$start_utc)) stop("Datetime input requires a dated study origin.", call. = FALSE)
  as.numeric(difftime(fr_parse_iso(x, name), origin$start_utc, units = "hours"))
}

fr_repair_headers <- function(d) {
  if (!ncol(d)) stop("The file contains no columns.", call. = FALSE)
  original <- names(d)
  repaired <- trimws(sub("^\ufeff", "", original))
  blank <- is.na(repaired) | !nzchar(repaired)
  repaired[blank] <- paste0("column_", which(blank))
  # make.unique reserves existing suffixes, so no column is merged or overwritten.
  repaired <- make.unique(repaired, sep = "__")
  names(d) <- repaired
  attr(d, "column_name_map") <- data.frame(position = seq_along(original),
    original = original, imported = repaired,
    changed = is.na(original) | original != repaired, stringsAsFactors = FALSE)
  d
}

fr_read_table <- function(path, filename, delimiter = "auto") {
  ext <- tolower(tools::file_ext(filename))
  # Examine the whole worksheet when determining types, avoiding late text
  # being silently coerced to missing. Native numeric Excel cells retain type.
  if (ext %in% c("xlsx", "xls")) d <- as.data.frame(readxl::read_excel(path, guess_max = 1048576, .name_repair = "minimal"))
  else if (ext %in% c("sav", "zsav")) d <- as.data.frame(haven::read_sav(path, .name_repair = "minimal"))
  else {
    delim <- switch(delimiter, auto = NULL, comma = ",", semicolon = ";", tab = "\t")
    imported <- suppressWarnings(readr::read_delim(path, delim = delim,
      col_types = readr::cols(.default = readr::col_character()), name_repair = "minimal",
      na = c("", "NA"), trim_ws = TRUE, show_col_types = FALSE, progress = FALSE))
    if (nrow(readr::problems(imported))) stop("Text file has inconsistent fields or quoting; correct the import problems before continuing.",call.=FALSE)
    d <- as.data.frame(imported)
  }
  fr_repair_headers(d)
}

fr_actigraphy <- function(times, states, epoch_h, end) {
  fr_scalar(epoch_h, "Epoch duration", 0, open_lo = TRUE)
  states <- fr_parse_numbers(states, name = "Actigraphy state", allow_missing = FALSE)
  if (!length(times) || any(!is.finite(times)) || any(!states %in% c(0,1)) || length(times) != length(states))
    stop("Actigraphy needs finite timestamps and one state (0 sleep / 1 wake) per epoch.", call. = FALSE)
  tol <- 1e-7
  if (abs(times[1]) > tol || abs(tail(times,1)+epoch_h-end) > tol)
    stop("Epoch-start actigraphy must cover the entire study, from hour 0 through the confirmed final epoch end.", call. = FALSE)
  if (length(times)>1L && any(abs(diff(times)-epoch_h)>tol))
    stop("Actigraphy has duplicate, unordered, missing or irregular epochs; resolve gaps before simulation.", call. = FALSE)
  # Snap only timestamp roundoff at known epoch boundaries.
  times <- seq(0, by = epoch_h, length.out = length(times))
  runs <- rle(states); last <- cumsum(runs$lengths); first <- c(1L,head(last,-1L)+1L)
  asleep <- which(runs$values == 0)
  x <- data.frame(onset_h = times[first[asleep]], offset_h = pmin(end,times[last[asleep]]+epoch_h))
  fr_intervals(x, end)
}

fr_time_token <- function(name) {
  s <- tolower(gsub("[[:space:]]+", "", name)); day <- NA_real_
  pat <- "(dag|day)[_.-]?([0-9]+)|(?:^|[_.-])(?:d|t)([0-9]+)(?=[_.-]|$)"
  m <- regmatches(s,regexec(pat,s,perl=TRUE))[[1]]
  if (length(m)) { day <- as.numeric(if (nzchar(m[3])) m[3] else m[4]); s <- sub(pat,"_",s,perl=TRUE) }
  get_clock <- function(h,m) if (h<24 && m<60) h+m/60 else NA_real_
  m <- regmatches(s,regexec("([0-9]{1,2})[:hu.]([0-9]{2})(?![0-9])",s,perl=TRUE))[[1]]
  if (length(m)) return(list(kind="clock",clock=get_clock(as.numeric(m[2]),as.numeric(m[3])),day=day))
  m <- regmatches(s,regexec("(?:^|[^0-9])([0-9]{4})(?:[^0-9]|$)",s,perl=TRUE))[[1]]
  if (length(m)) return(list(kind="clock",clock=get_clock(as.numeric(substr(m[2],1,2)),as.numeric(substr(m[2],3,4))),day=day))
  m <- regmatches(s,regexec("([0-9]{1,2})[hu](?![0-9])",s,perl=TRUE))[[1]]
  if (length(m)) return(list(kind="clock",clock=get_clock(as.numeric(m[2]),0),day=day))
  if (grepl("^[0-9]{1,2}$",s)) return(list(kind="clock",clock=get_clock(as.numeric(s),0),day=day))
  # Last resort: any number anywhere in the name (not just a strict
  # prefix_number shape), so names like "Lapses_10_SD" or "measurement2b"
  # still resolve instead of falling through to "unknown". The LAST number
  # found is used, matching the common prefix_N convention.
  nums <- regmatches(s,gregexpr("[0-9]+",s))[[1]]
  if (length(nums)) return(list(kind="index",clock=as.numeric(nums[length(nums)]),day=day))
  list(kind="unknown",clock=NA_real_,day=day)
}

fr_propose_times <- function(columns, t0_clock, step = 8) {
  fr_scalar(step,"Proposed occasion spacing",0,open_lo=TRUE)
  p <- lapply(columns,fr_time_token)
  if (!length(p)) return(setNames(numeric(),character()))
  kind <- vapply(p,`[[`,character(1),"kind"); cc <- vapply(p,`[[`,numeric(1),"clock")
  dd <- vapply(p,`[[`,numeric(1),"day"); ans <- rep(NA_real_,length(p))
  if (all(kind=="index") && all(is.na(dd))) ans <- (cc-1)*step
  else if (all(kind=="clock") && all(is.finite(cc))) {
    if (all(is.finite(dd)) && all(dd>=1)) ans <- (dd-1)*24+cc-t0_clock
    else if (all(is.na(dd))) {
      ans <- cc-t0_clock
      if (ans[1]<0) ans[1] <- ans[1]+24
      if (length(ans)>1L) for (i in 2:length(ans)) {
        if (cc[i]==cc[i-1L]) return(setNames(rep(NA_real_,length(ans)),columns))
        while (ans[i]<ans[i-1L]) ans[i] <- ans[i]+24
      }
    }
  }
  setNames(ans,columns)
}

# Bulk assignment uses the selected column order, never alphabetical sorting.
fr_bulk_times <- function(columns, mode = "regular", first = 0, step = 2,
                          t0_clock = 8, first_day = 1) {
  if (!length(columns)) stop("Select outcome columns first.",call.=FALSE)
  if (anyDuplicated(columns)) stop("Selected columns must be unique.",call.=FALSE)
  if (mode == "regular") {
    fr_scalar(first,"First measurement (hours after study start)",0)
    fr_scalar(step,"Time between measurements",0,open_lo=TRUE)
    ans <- first + (seq_along(columns)-1)*step
  } else if (mode == "clock") {
    fr_scalar(t0_clock,"Study-start clock",0,24,open_hi=TRUE)
    fr_scalar(first_day,"Day of first measurement",1)
    if (first_day != floor(first_day)) stop("Day of first measurement must be a whole number.",call.=FALSE)
    clocks <- vapply(columns,function(nm) {
      # This mode explicitly treats a number found anywhere in the column
      # name as a clock hour, regardless of prefix/suffix text (e.g.
      # "Lapses_10_SD" -> 10). fr_time_token() already tries clock-specific
      # patterns (H:MM, HHMM, NNh) before falling back to "any number", and
      # its own day-prefix handling is reused as-is (reject if a day marker
      # was found, since first_day above already covers multi-day numbering).
      token<-fr_time_token(nm)
      if(is.finite(token$clock) && is.na(token$day)) token$clock else NA_real_
    },numeric(1))
    if(any(!is.finite(clocks))) stop("Could not read a clock time from: ",
      paste(columns[!is.finite(clocks)],collapse=", "),
      ". Use regular spacing or the advanced day/time suggestion.",call.=FALSE)
    ans<-clocks-t0_clock+(first_day-1)*24
    # Place the first measurement on/after the chosen 24-hour study day;
    # later measurements use the next occurrence of their clock time.
    if(ans[1]<(first_day-1)*24) ans[1]<-ans[1]+24
    if(length(ans)>1L) for(i in 2:length(ans)) {
      if(ans[i]<=ans[i-1L]) ans[i]<-ans[i]+24*(floor((ans[i-1L]-ans[i])/24)+1)
    }
  } else stop("Choose regular spacing or clock times in column names.",call.=FALSE)
  if(any(!is.finite(ans))) stop("Measurement times are too large.",call.=FALSE)
  setNames(ans,columns)
}

fr_prepare_observations <- function(raw, settings, metric) {
  columns <- settings$columns
  if (!length(columns) || any(!columns %in% names(raw))) stop(metric, ": select outcome columns.",call.=FALSE)
  if (length(settings$times)!=length(columns) || any(!is.finite(settings$times)))
    stop(metric, ": use Fill all measurement times in the observations panel, or complete the remaining blank times in the table.",call.=FALSE)
  if (anyDuplicated(settings$times)) stop(metric, ": duplicate measurement times; aggregate within-subject replicates first.",call.=FALSE)
  if (settings$id == "__row__") subject <- as.character(seq_len(nrow(raw)))
  else {
    if (!settings$id %in% names(raw)) stop(metric, ": select a subject identifier.",call.=FALSE)
    subject <- trimws(as.character(raw[[settings$id]]))
  }
  group <- if (settings$group == "__none__") rep("All", nrow(raw)) else trimws(as.character(raw[[settings$group]]))
  if (length(group)!=nrow(raw)) stop(metric, ": missing subject/group labels must be resolved.",call.=FALSE)
  # Rows with a blank subject or group label are excluded rather than
  # aborting the whole import; missing category labels are common in
  # real datasets (e.g. a covariate that wasn't scored for everyone).
  keep <- !is.na(subject) & nzchar(subject) & !is.na(group) & nzchar(group)
  if (!all(keep)) { raw <- raw[keep, , drop=FALSE]; subject <- subject[keep]; group <- group[keep] }
  if (!nrow(raw)) stop(metric, ": no rows have both a subject and a group label.",call.=FALSE)
  if (anyDuplicated(data.frame(subject,group))) stop(metric, ": one row per subject within each group is required.",call.=FALSE)
  if (any(columns %in% c(settings$id,settings$group))) stop(metric, ": identifier/group columns cannot be outcome columns.",call.=FALSE)
  if (metric=="PVT" && !isTRUE(settings$pvt_units_confirmed))
    stop("Confirm that the PVT columns contain lapse counts from 10-minute tests.",call.=FALSE)
  if (metric=="PVT" && !settings$lapse_rule %in% c(">=500 ms",">500 ms")) stop("Specify the PVT lapse rule.",call.=FALSE)
  pieces <- lapply(seq_along(columns),function(j) {
    v <- fr_parse_numbers(raw[[columns[j]]],settings$decimal,paste(metric,columns[j]))
    data.frame(subject=subject,group=group,occasion=columns[j],time_h=settings$times[j],value_raw=v,value=v)
  })
  d <- do.call(rbind,pieces); rownames(d)<-NULL
  if (metric=="KSS" && any(d$value<1 | d$value>9,na.rm=TRUE)) stop("KSS observations must be between 1 and 9.",call.=FALSE)
  if (metric=="PVT" && any(d$value<0,na.rm=TRUE)) stop("PVT lapse counts cannot be negative.",call.=FALSE)
  if (metric=="VAS") {
    fr_scalar(settings$vas_min,"VAS minimum");fr_scalar(settings$vas_max,"VAS maximum",settings$vas_min,open_lo=TRUE)
    if (any(d$value<settings$vas_min | d$value>settings$vas_max,na.rm=TRUE))
      stop("VAS values outside the specified endpoints must be corrected; they are not silently clipped.",call.=FALSE)
    f <- (d$value-settings$vas_min)/(settings$vas_max-settings$vas_min)
    if (settings$direction=="reverse") f <- 1-f
    d$value <- if (identical(settings$target,"01")) f else 1+8*f
  }
  d$metric <- metric
  d
}

# Lin's concordance correlation coefficient (1989): rho_c = 2*cov(p,o) /
# (var(p)+var(o)+(mean(p)-mean(o))^2), using population (n-divisor) moments
# per the original definition. Combines precision (Pearson r) and accuracy
# (closeness to the 1:1 line) into one number in [-1,1]: 1 = perfect
# agreement, 0 = none, negative = systematic disagreement. This is the
# agreement statistic actually reported in sleep/chronobiology and clinical
# measurement-comparison literature -- unlike Nash-Sutcliffe efficiency (the
# statistic this replaces here), it's bounded and is not confounded by how
# much natural variance the particular sample happens to have.
fr_ccc <- function(p, o) {
  ok <- is.finite(p)&is.finite(o); p<-p[ok];o<-o[ok];n<-length(p)
  if (n<2) return(NA_real_)
  mp<-mean(p);mo<-mean(o)
  denom <- mean((p-mp)^2)+mean((o-mo)^2)+(mp-mo)^2
  if (denom<1e-12) return(NA_real_)
  2*mean((p-mp)*(o-mo))/denom
}

# Inter-observer CCC between two genuinely subject-varying series (e.g. two
# loaded outcome types, not a model prediction), following the GEE approach
# of Tsai (2017, Statist Med, doi:10.1002/sim.7207), Eq. (4)/(9), with the
# identity (D1) time-weight matrix -- the paper's own recommended default
# absent a specific reason to weight occasions unevenly. With no covariates,
# GEE's mean/variance/covariance estimating equations (Section 2.3) reduce
# exactly to simple sample moments regardless of the working-covariance
# choice, so no GEE-fitting package is needed here; only the pooling order
# (sum numerator/denominator components PER OCCASION, then ratio) is the
# part that must be followed precisely -- that's what distinguishes this
# from a naive flat CCC across all pooled (subject, occasion) pairs, and
# what stops a between-occasion mean shift from leaking into the estimate.
# NOTE: this is NOT the right tool for "model vs data" agreement -- a
# deterministic model prediction has zero between-subject variance at any
# fixed occasion, which forces Cov(model, data) = 0 by construction (see
# fr_ccc()/comp_data() for the model-vs-data statistics instead). Only use
# this between two series that both vary genuinely by subject.
fr_inter_ccc_gee <- function(subject1, time1, value1, subject2, time2, value2) {
  d1 <- data.frame(subject = as.character(subject1), time_h = time1, v1 = value1)
  d2 <- data.frame(subject = as.character(subject2), time_h = time2, v2 = value2)
  m <- merge(d1, d2, by = c("subject","time_h"))
  m <- m[is.finite(m$v1) & is.finite(m$v2), , drop = FALSE]
  empty <- list(ccc = NA_real_, n_pairs = 0L, n_times = 0L, n_subjects = 0L)
  if (!nrow(m)) return(empty)
  num <- 0; den <- 0; n_times <- 0L
  for (tk in sort(unique(m$time_h))) {
    z <- m[m$time_h == tk, ]
    if (nrow(z) < 2) next
    mu1 <- mean(z$v1); mu2 <- mean(z$v2)
    var1 <- mean((z$v1-mu1)^2); var2 <- mean((z$v2-mu2)^2)
    cov12 <- mean((z$v1-mu1)*(z$v2-mu2))
    num <- num + 2*cov12
    den <- den + var1 + var2 + (mu1-mu2)^2
    n_times <- n_times + 1L
  }
  ccc <- if (n_times >= 1 && den > 1e-12) num/den else NA_real_
  list(ccc = ccc, n_pairs = nrow(m), n_times = n_times, n_subjects = length(unique(m$subject)))
}

fr_fit_statistics <- function(p, o) {
  ok <- is.finite(p)&is.finite(o); p<-p[ok];o<-o[ok];n<-length(p)
  # Pearson_r only measures whether pred and obs move together (shape/rank);
  # it is blind to a systematic offset or scale mismatch, so a high r can
  # coexist with a poor fit. CCC is the statistic that should read as
  # "goodness of fit" -- see fr_ccc() above.
  c(pairs=n, RMSE=if(n) sqrt(mean((p-o)^2)) else NA_real_,
    MAE=if(n) mean(abs(p-o)) else NA_real_, bias=if(n) mean(p-o) else NA_real_,
    Pearson_r=if(n>=2 && sd(p)>1e-12 && sd(o)>1e-12) cor(p,o) else NA_real_,
    CCC=fr_ccc(p,o))
}

fr_compare <- function(pred, observations, selected_groups, schedule, end, metric, selected_subjects = NULL) {
  d <- observations; d$exclusion <- ""
  d$exclusion[is.na(d$value)] <- "Missing observation"
  d$exclusion[d$time_h<0 | d$time_h>end] <- "Outside study coverage"
  d$exclusion[fr_inside(d$time_h,schedule)] <- "During sleep"
  d$exclusion[!d$group %in% selected_groups] <- "Group not selected"
  if (!is.null(selected_subjects)) d$exclusion[!as.character(d$subject) %in% selected_subjects] <- "Participant not selected"
  good <- d[d$exclusion=="",,drop=FALSE]
  cells <- data.frame(group=character(),occasion=character(),time_h=numeric(),n=integer(),
    observed=numeric(),SD=numeric(),SEM=numeric(),pred_raw=numeric(),pred_bounded=numeric(),
    residual_raw=numeric(),residual_bounded=numeric())
  if (nrow(good)) {
    keys <- interaction(good$group,good$occasion,drop=TRUE,lex.order=TRUE)
    cells <- do.call(rbind,lapply(split(good,keys),function(z) {
      idx <- match(z$time_h[1],pred$time_h)
      if (is.na(idx)) stop("An eligible observation time is missing from the solver mesh.",call.=FALSE)
      # VAS01: a VAS/other-scale import rescaled to 0-1 (rather than to the
      # 1-9 KSS-equivalent default) compares against the model's own composite
      # risk index (Index_latent/Index), which is always computed regardless
      # of what empirical data is loaded -- unlike Index_compare, which needs
      # a KSS+PVT-adaptive index_settings entry that a standalone VAS-only
      # import wouldn't have.
      col_raw <- switch(metric, PVT = "PVT_raw", Index = "Index_compare_raw", VAS01 = "Index_latent", "KSS_raw")
      col_bounded <- switch(metric, PVT = "PVT", Index = "Index_compare", VAS01 = "Index", "KSS")
      p <- pred[[col_raw]][idx];b<-pred[[col_bounded]][idx]
      o <- mean(z$value); s<-sd(z$value);n<-nrow(z)
      data.frame(group=z$group[1],occasion=z$occasion[1],time_h=z$time_h[1],n=n,observed=o,SD=s,SEM=s/sqrt(n),
        pred_raw=p,pred_bounded=b,residual_raw=p-o,residual_bounded=b-o)
    }))
    cells <- cells[order(cells$group,cells$time_h),,drop=FALSE];rownames(cells)<-NULL
  }
  scopes <- c("Pooled",sort(unique(cells$group)))
  stats <- do.call(rbind,lapply(seq_along(scopes),function(i) {
    z <- if (i==1L) cells else cells[cells$group==scopes[i],,drop=FALSE]
    do.call(rbind,lapply(c("raw","bounded"),function(scale) {
      values<-fr_fit_statistics(z[[paste0("pred_",scale)]],z$observed)
      data.frame(scope=if(i==1L)"Pooled (equal cell weights)" else paste0("Group: ",scopes[i]),
                 prediction_scale=scale,as.list(values),check.names=FALSE)
    }))
  }))
  list(cells=cells,statistics=stats,observations=d,
       exclusions=d[d$exclusion!="",,drop=FALSE],selected_groups=selected_groups,selected_subjects=selected_subjects)
}

fr_run <- function(inputs) {
  end <- fr_scalar(inputs$end,"Study duration",0,open_lo=TRUE)
  schedule <- fr_intervals(inputs$schedule,end)
  duty <- fr_intervals(inputs$duty,end,"Duty schedule")
  if (inputs$exposure$scope=="duty" && !nrow(duty)) stop("Enter duty intervals for duty-scoped exposure.",call.=FALSE)
  I <- fr_index_settings(inputs$index)
  extra <- unlist(lapply(inputs$empirical,function(e) e$observations$time_h),use.names=FALSE)
  extra <- c(extra,unlist(duty,use.names=FALSE))
  arms <- lapply(c("PVT","KSS"),function(m) {
    P <- fr_parameters(m,inputs$version,inputs$overrides[[m]])
    fr_simulate_arm(P,schedule,end,inputs$origin$t0_clock,inputs$spacing,inputs$baseline,inputs$initial[[m]],extra)
  }); names(arms)<-c("PVT","KSS")
  a<-arms$PVT$predictions;b<-arms$KSS$predictions
  stopifnot(identical(a$time_h,b$time_h))
  pred<-a[,c("time_h","clock_h","wake","c")]
  pred$PVT_raw<-a$p;pred$KSS_raw<-b$p;pred$PVT<-pmax(a$p,0);pred$KSS<-pmin(pmax(b$p,1),9)
  pred$q_PVT<-a$q;pred$q_KSS<-b$q
  pred$u_PVT<-a$u;pred$u_KSS<-b$u;pred$h_PVT<-a$h;pred$h_KSS<-b$h
  pred$kappa_PVT<-a$kappa;pred$kappa_KSS<-b$kappa
  pred$inertia_PVT<-fr_parameters("PVT",inputs$version,inputs$overrides$PVT)$xi*a$h
  pred$inertia_KSS<-fr_parameters("KSS",inputs$version,inputs$overrides$KSS)$xi*b$h
  pred$Index_latent<-fr_index(pred$PVT_raw,pred$KSS_raw,I)
  pred$Index<-ifelse(pred$wake,pred$Index_latent,NA_real_)
  pred$on_duty<-fr_inside(pred$time_h,duty)
  # The composite comparison uses whichever weighting is actually testable: the
  # configured weights when both KSS and PVT observations are available, or
  # 100% on the single available term when only one is -- so the predicted
  # series it is scored against matches what was actually observed, rather
  # than always scoring against the globally configured (e.g. mixed) Index.
  if (!is.null(inputs$empirical$Index$index_settings)) {
    Ic <- inputs$empirical$Index$index_settings
    pred$Index_compare_raw <- fr_index(pred$PVT_raw, pred$KSS_raw, Ic)
    pred$Index_compare <- ifelse(pred$wake, pred$Index_compare_raw, NA_real_)
  }
  exposure<-fr_exposure(pred,schedule,duty,inputs$exposure,I)
  comparisons<-lapply(names(inputs$empirical),function(m) {
    e<-inputs$empirical[[m]]
    fr_compare(pred,e$observations,e$selected_groups,schedule,end,e$compare_as %||% m,e$selected_subjects)
  });names(comparisons)<-names(inputs$empirical)
  P<-setNames(lapply(c("PVT","KSS"),function(m)fr_parameters(m,inputs$version,inputs$overrides[[m]])),c("PVT","KSS"))
  result<-list(predictions=pred,exposure=exposure,comparisons=comparisons,
    parameters=P,initialization=lapply(arms,`[[`,"initialization"),
    range_counts=data.frame(outcome=c("PVT","KSS"),
      waking_predictions_outside_scale=c(sum(pred$wake & pred$PVT_raw<0),sum(pred$wake & (pred$KSS_raw<1|pred$KSS_raw>9)))))
  pkgs<-c("deSolve","shiny","shinydashboard","DT","plotly","readr","readxl","haven","lubridate","digest")
  list(inputs=inputs,results=result,metadata=list(app_version=FR_VERSION,created_utc=format(Sys.time(),tz="UTC",usetz=TRUE),
    R=R.version.string,packages=setNames(vapply(pkgs,function(p) if(requireNamespace(p,quietly=TRUE))as.character(utils::packageVersion(p)) else NA_character_,character(1)),pkgs),
    references=FR_REFERENCES,exposure_method="Piecewise-linear positive-part integration with threshold and clipping crossings; adaptive ODE states; a threshold-exposure adaptation, not a calibrated accident probability"))
}

# Searches the checked parameters (bounded the same way as their manual
# sliders) to maximize CCC (Lin's concordance correlation coefficient, see
# fr_ccc()) against the currently loaded, currently selected observations for
# one metric. Plain function -- no Shiny dependency -- so it can be called
# repeatedly by an optimizer outside the app's normal debounced reactive
# pipeline. CCC in [-1,1] also gives optim() a bounded objective, unlike the
# previous NSE-based objective which was unbounded below.
fr_suggest_calibration <- function(base_inputs, arm, comparison_metric, param_names,
                                    scale = "raw", maxit = 150) {
  stopifnot(arm %in% c("PVT","KSS"))
  cur <- fr_parameters(arm, base_inputs$version, base_inputs$overrides[[arm]])
  param_names <- intersect(param_names, names(cur))
  if (!length(param_names)) stop("No adjustable parameters selected.", call. = FALSE)
  x0 <- unlist(cur[param_names])
  bounds <- vapply(param_names, function(nm) fr_slider_range(nm, cur[[nm]]), numeric(2))
  lower <- bounds[1,]; upper <- bounds[2,]

  eval_ccc <- function(par_vec) {
    ov <- as.list(par_vec); names(ov) <- param_names
    inputs2 <- base_inputs
    inputs2$overrides[[arm]] <- utils::modifyList(inputs2$overrides[[arm]] %||% list(), ov)
    tryCatch({
      cells <- fr_run(inputs2)$results$comparisons[[comparison_metric]]$cells
      if (is.null(cells) || nrow(cells) < 2) return(NA_real_)
      fr_fit_statistics(cells[[paste0("pred_", scale)]], cells$observed)[["CCC"]]
    }, error = function(e) NA_real_)
  }
  objective <- function(par_vec) { ccc <- eval_ccc(par_vec); if (!is.finite(ccc)) 1e6 else -ccc }

  ccc_before <- eval_ccc(x0)
  fit <- optim(x0, objective, method = "L-BFGS-B", lower = lower, upper = upper,
               control = list(maxit = maxit, parscale = pmax(abs(x0), 1e-3)))
  ccc_after <- -fit$value
  at_bound <- (abs(fit$par - lower) < 1e-6) | (abs(fit$par - upper) < 1e-6)
  # Strip names here so callers never round-trip a named atomic vector into a
  # Shiny client message (e.g. updateSliderInput(value=...)) -- jsonlite's
  # asJSON(keep_vec_names=TRUE) turns a named length-1 vector into a JSON
  # object instead of a plain value, which malforms that client message.
  list(arm = arm, param_names = param_names, before = unname(x0), after = unname(fit$par),
       lower = unname(lower), upper = unname(upper), at_bound = unname(at_bound),
       ccc_before = ccc_before, ccc_after = ccc_after,
       convergence = fit$convergence, message = fit$message)
}

fr_code_snapshot <- function() {
  setNames(vapply(FR_CORE_FUNCTIONS,function(nm)paste(deparse(get(nm,envir=environment(fr_run)),width.cutoff=500),collapse="\n"),character(1)),FR_CORE_FUNCTIONS)
}

fr_capture_run <- function(inputs) {
  b <- fr_run(inputs)
  b$code <- fr_code_snapshot()
  b$metadata$code_sha256 <- digest::digest(b$code,algo="sha256")
  b
}

fr_export_script <- function(bundle) {
  dump <- function(x) paste(capture.output(dput(x,control="all")),collapse="\n")
  definitions <- paste(paste(names(bundle$code),"<-",bundle$code),collapse="\n\n")
  paste0("# FaRiMaTo reproducible analysis — ",bundle$metadata$app_version,"\n",
    "# Generated from the same immutable inputs and processed observations as the displayed run.\n",
    "# Package versions are checked below; no packages are installed automatically.\n",
    "FR_VERSION <- ",dump(bundle$metadata$app_version),"\nFR_REFERENCES <- ",dump(bundle$metadata$references),
    "\nFR_CORE_FUNCTIONS <- ",dump(names(bundle$code)),"\n\n",definitions,
    "\n\nrecorded_packages <- ",dump(bundle$metadata$packages),
    "\nfor (p in names(recorded_packages)) {\n",
    "  if (!is.na(recorded_packages[[p]]) && (!requireNamespace(p, quietly=TRUE) || as.character(utils::packageVersion(p)) != recorded_packages[[p]]))\n",
    "    stop(paste('Reproduction requires recorded package',p,recorded_packages[[p]]))\n}\n",
    "\ninputs <- ",dump(bundle$inputs),"\nexpected <- ",dump(bundle$results),
    "\nreproduced <- fr_run(inputs)\n",
    "check <- all.equal(reproduced$results, expected, tolerance=1e-7)\n",
    "if (!isTRUE(check)) stop(paste('Reproduction mismatch:',paste(check,collapse='; ')))\n",
    "saveRDS(reproduced, 'farimato_reproduced.rds')\n",
    "write.csv(reproduced$results$predictions,'farimato_reproduced_predictions.csv',row.names=FALSE)\n",
    "message('Verified: all prediction, exposure, initialization and comparison outputs reproduced.')\n")
}

# Captured before UI helpers, so standalone exports contain every computational
# dependency without importing a Shiny application or relying on a source path.
FR_CORE_FUNCTIONS <- c("fr_scalar", "mccauley_params", "fr_parameters", "fr_empty_intervals",
  "fr_intervals", "fr_inside", "fr_rhs", "fr_advance", "fr_endpoint", "fr_baseline",
  "fr_initial_state", "fr_simulate_arm", "fr_index_settings", "fr_index", "fr_positive_segment",
  "fr_exposure", "fr_parse_numbers", "fr_parse_times", "fr_parse_iso", "fr_origin", "fr_elapsed",
  "fr_repair_headers", "fr_read_table", "fr_actigraphy", "fr_time_token", "fr_propose_times", "fr_bulk_times", "fr_prepare_observations",
  "fr_fit_statistics", "fr_compare", "fr_run", "fr_code_snapshot", "fr_capture_run", "fr_export_script")

# Elapsed-hour tick positions relabelled as clock-of-day text (HH:MM), so plot
# x-axes keep a straight monotonic timeline but read as clock times. A fixed
# `step` (hours) overrides the default "about n round ticks" auto-spacing.
fr_clock_ticks <- function(end, t0_clock, step = NULL, n = 8) {
  if (!is.finite(end) || end <= 0) end <- 1
  if (!is.finite(t0_clock)) t0_clock <- 0
  vals <- if (!is.null(step) && is.finite(step) && step > 0) seq(0, end, by = step)
          else pretty(c(0, end), n = n)
  vals <- vals[vals >= 0 & vals <= end]
  if (!length(vals)) vals <- 0
  clock <- (t0_clock + vals) %% 24
  hh <- floor(clock + 1e-9); mm <- round((clock - hh) * 60)
  carry <- mm >= 60; mm[carry] <- 0; hh[carry] <- (hh[carry] + 1) %% 24
  list(tickvals = vals, ticktext = sprintf("%02d:%02d", hh, mm))
}

# Reads the "Time axis resolution" selectInput ("" = auto) into a numeric step.
fr_clock_step_from_input <- function(input) {
  v <- input$clock_tick_hours
  if (is.null(v) || !nzchar(v)) NULL else as.numeric(v)
}

# Shared clock-labelled x-axis for a predictions data frame carrying the run's
# bundle (t0_clock) as an attribute, as produced by the sim_results reactive.
fr_xaxis_clock <- function(res, step = NULL, extra = list()) {
  t0 <- attr(res, "bundle")$inputs$origin$t0_clock
  end <- suppressWarnings(max(res$time_h, na.rm = TRUE))
  ticks <- fr_clock_ticks(end, t0, step)
  utils::modifyList(list(title = "Clock time", zeroline = FALSE,
    tickvals = ticks$tickvals, ticktext = ticks$ticktext), extra)
}

hex_to_rgba <- function(hex, alpha = 0.25) {
  r <- strtoi(substr(hex, 2L, 3L), 16L)
  g <- strtoi(substr(hex, 4L, 5L), 16L)
  b <- strtoi(substr(hex, 6L, 7L), 16L)
  sprintf("rgba(%d,%d,%d,%.2f)", r, g, b, alpha)
}

# Helper: build a distinct colour per group (Set2 pastel, House Rule 9)
group_colours <- function(groups) {
  if(!length(groups)) return(setNames(character(),character()))
  cols<-if(length(groups)==1) "#e74c3c" else grDevices::colorRampPalette(RColorBrewer::brewer.pal(8,"Set2"))(length(groups))
  setNames(cols,groups)
}
add_outcome_overlay <- function(p, emp_df, value_col, yaxis = "y",
                                show_indiv = FALSE, show_mean = TRUE,
                                label_prefix = "Observed",
                                line_dash = "solid", marker_symbol = "circle",
                                y_clip = NULL, color_override = NULL) {
  if (is.null(emp_df) || nrow(emp_df) == 0 || !(show_indiv || show_mean)) return(p)

  groups   <- sort(unique(emp_df$group))
  n_groups <- length(groups)
  # When several outcome types (KSS/VAS/PVT) are overlaid on the same plot,
  # color_override fixes one colour per outcome so they're distinguishable at
  # a glance instead of all landing on the same first group colour; groups
  # within a single outcome still get distinct colours when no override is set.
  grp_cols <- if (!is.null(color_override)) setNames(rep(color_override, n_groups), groups) else group_colours(groups)

  if (show_indiv) {
    for (grp in groups) {
      d_g <- emp_df[emp_df$group == grp, ]
      x_v <- y_v <- c()
      for (sid in unique(d_g$subject)) {
        d_s <- d_g[d_g$subject == sid, ]
        d_s <- d_s[order(d_s$time_h), ]
        x_v <- c(x_v, d_s$time_h, NA_real_)
        y_v <- c(y_v, d_s[[value_col]], NA_real_)
      }
      lbl <- if (n_groups <= 1) paste0(label_prefix, " (individual)")
             else paste0(label_prefix, ": ", grp, " (individual)")
      p <- add_trace(p, inherit=FALSE, x = x_v, y = y_v, type = "scatter", mode = "lines",
                     yaxis = yaxis,
                     line = list(color = hex_to_rgba(unname(grp_cols[[grp]]), 0.30),
                                 width = 0.9, dash = line_dash),
                     name = lbl, legendgroup = paste0("indiv_", label_prefix, "_", grp),
                     showlegend = TRUE, hoverinfo = "none")
    }
  }

  if (show_mean) {
    for (grp in groups) {
      d_grp <- emp_df[emp_df$group == grp, ] %>%
        dplyr::group_by(time_h) %>%
        dplyr::summarise(m = mean(.data[[value_col]], na.rm = TRUE),
                         s = sd(.data[[value_col]],   na.rm = TRUE),
                         n = sum(is.finite(.data[[value_col]])), .groups = "drop") %>%
        dplyr::arrange(time_h)

      col <- unname(grp_cols[[grp]])
      lbl <- if (n_groups <= 1) paste0(label_prefix, " (mean +/- SD)")
             else paste0(label_prefix, ": ", grp, " (mean +/- SD)")

      upper <- d_grp$m + d_grp$s
      lower <- d_grp$m - d_grp$s
      if (FALSE) {
        upper <- pmin(upper, y_clip[2])
        lower <- pmax(lower, y_clip[1])
      }

      # Only draw SD ribbons where SD is defined; never join across missing SD.
      valid<-which(is.finite(upper)&is.finite(lower))
      if(length(valid))for(idx in split(valid,cumsum(c(TRUE,diff(valid)>1)))) {
        if(length(idx)<2)next
        p<-add_trace(p,x=c(d_grp$time_h[idx],rev(d_grp$time_h[idx])),
          y=c(upper[idx],rev(lower[idx])),inherit=FALSE,type="scatter",mode="none",
          fill="toself",fillcolor=hex_to_rgba(col,.15),yaxis=yaxis,
          legendgroup=paste0("mean_",label_prefix,"_",grp),showlegend=FALSE,hoverinfo="none")
      }

      p <- add_trace(p, inherit=FALSE, data = d_grp, x = ~time_h, y = ~m,
                     type = "scatter", mode = "lines+markers", yaxis = yaxis,
                     legendgroup = paste0("mean_", label_prefix, "_", grp),
                     name = lbl,
                     line   = list(color = col, width = 2, dash = line_dash),
                     marker = list(color = col, size = 7, symbol = marker_symbol,
                                   line = list(color = col, width = 2)),
                     customdata  = ~n,
                     hovertemplate = paste0("<b>t = %{x:.1f} h</b><br>",
                                            label_prefix, " = %{y:.2f}<br>",
                                            "n = %{customdata}<extra></extra>"))
    }
  }
  p
}

make_comparison_outputs <- function(input, output, session, prefix, label,
                                    emp_reactive, sim_results_fn,
                                    value_col, model_col, y_range,
                                    err_type_input = "err_type") {

  # ---- Group selector -------------------------------------------------
  # Fit statistics and every plot below are computed on the SELECTED groups
  # only, so deselecting a group genuinely re-estimates RMSE / MAE / r rather
  # than merely hiding a trace.
  available_groups<-reactive({
    emp<-emp_reactive();if(is.null(emp))character() else sort(unique(as.character(emp$group)))
  })
  group_choices<-reactiveVal(character())
  observe({g<-tryCatch(available_groups(),error=function(e)NULL);if(is.null(g))return();if(!identical(g,isolate(group_choices())))group_choices(g)})
  output[[paste0(prefix, "_ui_groups")]] <- renderUI({
    grps<-group_choices();req(length(grps))
    if (length(grps) <= 1) return(NULL)   # nothing to choose from
    pickerInput(paste0(prefix, "_groups_selected"),
                label    = "Groups included in plots and fit statistics:",
                choices  = grps, selected = grps, multiple = TRUE,
                options  = list(`actions-box` = TRUE,
                                `selected-text-format` = "count > 3",
                                `count-selected-text`  = "{0} of {1} groups"))
  })

  # ---- Participant selector -------------------------------------------
  # Same mechanism as the group selector, at the individual-subject level.
  # Defaults to every participant included; deselecting one re-estimates
  # RMSE / MAE / r and drops that participant from the overlay traces too.
  available_subjects<-reactive({
    emp<-emp_reactive();if(is.null(emp))character() else sort(unique(as.character(emp$subject)))
  })
  subject_choices<-reactiveVal(character())
  observe({s<-tryCatch(available_subjects(),error=function(e)NULL);if(is.null(s))return();if(!identical(s,isolate(subject_choices())))subject_choices(s)})
  output[[paste0(prefix, "_ui_subjects")]] <- renderUI({
    subs<-subject_choices();req(length(subs))
    if (length(subs) <= 1) return(NULL)   # nothing to choose from
    pickerInput(paste0(prefix, "_subjects_selected"),
                label    = "Participants included in plots and fit statistics:",
                choices  = subs, selected = subs, multiple = TRUE,
                options  = list(`actions-box` = TRUE, `live-search` = TRUE,
                                `selected-text-format` = "count > 3",
                                `count-selected-text`  = "{0} of {1} participants"))
  })

  comp_data <- reactive({
    res<-sim_results_fn();req(!is.null(res))
    metric<-switch(prefix,kss="KSS",vas="VAS",pvtlap="PVT",index="Index")
    z<-attr(res,"bundle")$results$comparisons[[metric]]
    req(!is.null(z))
    d<-z$cells;scale<-input$prediction_scale
    data.frame(time_h=d$time_h,group=d$group,obs_mean=d$observed,obs_sd=d$SD,
      obs_sem=d$SEM,obs_n=d$n,pred=d[[paste0("pred_",scale)]],
      residual=d[[paste0("residual_",scale)]],abs_error=abs(d[[paste0("residual_",scale)]]))
  })
  # Individual-level pairing: every raw (subject, occasion) observation kept
  # by fr_compare's exclusion filter, matched to the SAME cell-level model
  # prediction its group/occasion already has (the model doesn't vary by
  # subject). Used for the individual-level CCC, so between-subject noise
  # that the group-mean cells average away counts against agreement here.
  comp_data_indiv <- reactive({
    res<-sim_results_fn();req(!is.null(res))
    metric<-switch(prefix,kss="KSS",vas="VAS",pvtlap="PVT",index="Index")
    z<-attr(res,"bundle")$results$comparisons[[metric]]
    req(!is.null(z))
    good<-z$observations[z$observations$exclusion=="",,drop=FALSE]
    req(nrow(good)>0)
    cellp<-z$cells[,c("group","occasion","pred_raw","pred_bounded")]
    d<-merge(good,cellp,by=c("group","occasion"))
    scale<-input$prediction_scale
    data.frame(time_h=d$time_h,group=d$group,subject=d$subject,
      obs=d$value,pred=d[[paste0("pred_",scale)]])
  })
  output[[paste0(prefix,"_exclusions")]]<-renderText({
    res<-sim_results_fn();req(!is.null(res))
    metric<-switch(prefix,kss="KSS",vas="VAS",pvtlap="PVT",index="Index")
    z<-attr(res,"bundle")$results$comparisons[[metric]];req(!is.null(z))
    counts<-table(z$exclusions$exclusion)
    paste("Comparison uses",nrow(z$cells),"group/occasion means with equal weights.",
      if(length(counts))paste("Excluded:",paste(names(counts),counts,collapse="; ")) else "No observations excluded.",
      if(metric=="VAS")"Rescaling from the raw scale is exploratory, not a validated conversion." else "")
  })

  fit_stat <- function(fn) {
    tryCatch({
      cd <- comp_data(); req(!is.null(cd) && nrow(cd) >= 1)
      d <- cd[!is.na(cd$pred) & !is.na(cd$obs_mean), ]
      if (nrow(d) < 1) return(NA_real_)
      fn(d$pred, d$obs_mean)
    }, error = function(e) NA_real_)
  }
  fit_stat_indiv <- function(fn) {
    tryCatch({
      cd <- comp_data_indiv(); req(!is.null(cd) && nrow(cd) >= 1)
      d <- cd[!is.na(cd$pred) & !is.na(cd$obs), ]
      if (nrow(d) < 1) return(NA_real_)
      fn(d$pred, d$obs)
    }, error = function(e) NA_real_)
  }

  span <- diff(y_range)

  # Subtitle suffix naming the groups the statistic was computed on
  grp_suffix <- reactive({
    cd <- tryCatch(comp_data(), error = function(e) NULL)
    if (is.null(cd) || nrow(cd) == 0) return("")
    g <- sort(unique(as.character(cd$group)))
    if (length(g) == 1 && identical(g, "All")) return("")
    if (length(g) <= 2) paste0(" | ", paste(g, collapse = " + "))
    else paste0(" | ", length(g), " groups")
  })

  output[[paste0(prefix, "_vb_rmse")]] <- renderValueBox({
    val <- fit_stat(function(p, o) sqrt(mean((p - o)^2)))
    col <- "blue"
    valueBox(if (is.na(val)) "--" else sprintf("%.2f", val),
             paste0("RMSE (", label, ")", grp_suffix()),
             icon = icon("bullseye"), color = col)
  })
  output[[paste0(prefix, "_vb_mae")]] <- renderValueBox({
    val <- fit_stat(function(p, o) mean(abs(p - o)))
    col <- "blue"
    valueBox(if (is.na(val)) "--" else sprintf("%.2f", val),
             paste0("MAE (", label, ")", grp_suffix()),
             icon = icon("ruler-combined"), color = col)
  })
  output[[paste0(prefix, "_vb_r")]] <- renderValueBox({
    val <- fit_stat(function(p, o) if (length(p)<2 || sd(p) < 1e-9 || sd(o) < 1e-9) NA_real_ else cor(p, o))
    col <- "blue"
    valueBox(if (is.na(val)) "--" else sprintf("%.3f", val),
             paste0("Pearson r", grp_suffix()),
             icon = icon("chart-line"), color = col)
  })
  # Lin's concordance correlation coefficient (fr_ccc): unlike Pearson r this
  # is penalized by a systematic offset/scale mismatch, so it can go negative
  # -- that is the statistic to read as "does predicted actually agree with
  # observed". Reported at two levels of aggregation (see helpText below each
  # box): group-level uses the same group/occasion-mean cells as RMSE/MAE/r
  # above; individual-level uses every raw kept observation instead, so
  # between-subject variability that the group means average away counts
  # against it -- typically lower, and a large gap between the two indicates
  # the model tracks the average trend but not individual differences.
  output[[paste0(prefix, "_vb_ccc")]] <- renderValueBox({
    val <- fit_stat(fr_ccc)
    col <- if (!is.na(val) && val < 0) "red" else "blue"
    valueBox(if (is.na(val)) "--" else sprintf("%.3f", val),
             paste0("CCC (group)", grp_suffix()),
             icon = icon("bullseye"), color = col)
  })
  output[[paste0(prefix, "_vb_ccc_ind")]] <- renderValueBox({
    val <- fit_stat_indiv(fr_ccc)
    col <- if (!is.na(val) && val < 0) "red" else "blue"
    valueBox(if (is.na(val)) "--" else sprintf("%.3f", val),
             paste0("CCC (individual)", grp_suffix()),
             icon = icon("users"), color = col)
  })
  # Per-subject CCC (each subject's own repeated occasions vs the model,
  # then averaged across subjects) -- built on the SAME filtered comp_data_indiv()
  # used above, so it respects the current group/subject selection exactly
  # like every other box here. A subject needs >=2 kept occasions with
  # nonzero variance in both series to contribute; fr_ccc() already returns
  # NA otherwise, so those subjects are dropped from the mean, not zeroed.
  # NOTE (see helpText below): this simple mean-of-independent-CCCs is what
  # was asked for, but Carrasco/King/Chinchilli-type variance-components or
  # GEE estimators (e.g. PMID 28105775) are the more rigorous way to pool
  # repeated-measures agreement across subjects, since they model the
  # within-subject correlation structure instead of averaging unweighted,
  # independently-estimated per-subject CCCs.
  ccc_by_subject <- reactive({
    cd <- tryCatch(comp_data_indiv(), error = function(e) NULL)
    empty <- list(mean = NA_real_, n_used = 0L, n_total = 0L)
    if (is.null(cd) || !nrow(cd)) return(empty)
    d <- cd[!is.na(cd$pred) & !is.na(cd$obs), ]
    if (!nrow(d)) return(empty)
    per_subj <- vapply(split(d, d$subject), function(s) fr_ccc(s$pred, s$obs), numeric(1))
    ok <- is.finite(per_subj)
    list(mean = if (any(ok)) mean(per_subj[ok]) else NA_real_,
         n_used = sum(ok), n_total = length(per_subj))
  })
  output[[paste0(prefix, "_vb_ccc_subj")]] <- renderValueBox({
    r <- ccc_by_subject(); val <- r$mean
    col <- if (!is.na(val) && val < 0) "red" else "blue"
    valueBox(if (is.na(val)) "--" else sprintf("%.3f", val),
             paste0("CCC (per-subject mean, n=", r$n_used, "/", r$n_total, ")", grp_suffix()),
             icon = icon("user"), color = col)
  })

  output[[paste0(prefix, "_plt_ts")]] <- renderPlotly({
    res <- sim_results_fn(); req(!is.null(res))
    cd  <- tryCatch(comp_data(), error = function(e) NULL)
    req(!is.null(cd) && nrow(cd) > 0)
    # A VAS import rescaled to 0-1 is compared against the composite risk
    # index (Risk, aliased from Index in sim_results), not the KSS model line.
    mcol <- if (identical(prefix,"vas") && identical(input$vas_target,"01")) "Risk" else model_col
    wake <- res; wake[[mcol]][!wake$wake] <- NA_real_

    groups   <- sort(unique(cd$group))
    grp_cols <- group_colours(groups)
    err_col  <- if (isTRUE(input[[err_type_input]] == "sem")) "obs_sem" else "obs_sd"

    p <- plot_ly(wake, x = ~time_h, y = wake[[mcol]],
                type = "scatter", mode = "lines",
                name = paste0(label, " model"),
                line = list(color = "#2980b9", width = 2.5))

    for (grp in groups) {
      d_g <- cd[cd$group == grp, ]
      col <- unname(grp_cols[[grp]])
      err <- d_g[[err_col]]
      lbl <- if (length(groups) <= 1) "Observed mean" else paste0("Observed: ", grp)
      p <- add_trace(p, data = d_g, x = ~time_h, y = ~obs_mean,
                     type = "scatter", mode = "markers+lines",
                     error_y = list(type = "data", array = err, color = col,
                                    thickness = 2, width = 8),
                     name = lbl, marker = list(color = col, size = 10),
                     line = list(color = col, width = 1.5, dash = "dot"),
                     text = ~paste0("t = ", round(time_h, 1), " h<br>Obs = ", round(obs_mean, 2),
                                    " +/- ", round(err, 2), " (n = ", obs_n, ")<br>",
                                    "Model = ", round(pred, 2), "<br>Residual = ", round(residual, 2)),
                     hoverinfo = "text")
    }
    p %>% layout(xaxis = fr_xaxis_clock(res,fr_clock_step_from_input(input)),
                yaxis = list(title = label, autorange = TRUE, gridcolor = "#eee"),
                shapes = sleep_shapes(res),
                legend = list(orientation = "h", y = -0.28),
                margin = list(t = 5, r = 10),
                plot_bgcolor = "#fefefe", paper_bgcolor = "#fefefe")
  })

  output[[paste0(prefix, "_plt_scatter")]] <- renderPlotly({
    cd <- tryCatch(comp_data(), error = function(e) NULL)
    req(!is.null(cd) && nrow(cd) >= 1)
    groups   <- sort(unique(cd$group))
    grp_cols <- group_colours(groups)
    err_col  <- if (isTRUE(input[[err_type_input]] == "sem")) "obs_sem" else "obs_sd"

    vals  <- c(cd$pred, cd$obs_mean)
    ax_lo <- floor(min(vals, na.rm = TRUE) - 0.05 * max(diff(range(vals)),1))
    ax_hi <- ceiling(max(vals, na.rm = TRUE) + 0.05 * max(diff(range(vals)),1))
    ax_range <- c(ax_lo, ax_hi)

    p <- plot_ly() %>%
      add_lines(x = ax_range, y = ax_range, name = "y = x (perfect fit)",
               line = list(color = "#999", dash = "dash", width = 1.5))

    for (grp in groups) {
      d_g <- cd[cd$group == grp & !is.na(cd$pred), ]
      col <- unname(grp_cols[[grp]])
      err <- d_g[[err_col]]
      lbl <- if (length(groups) <= 1) "Occasions" else grp
      p <- add_trace(p, data = d_g, x = ~pred, y = ~obs_mean,
                     type = "scatter", mode = "markers",
                     error_y = list(type = "data", array = err, color = col, thickness = 2, width = 5),
                     name = lbl, marker = list(color = col, size = 10,
                                               line = list(color = "white", width = 1)),
                     text = ~paste0("t = ", round(time_h, 1), " h<br>Model = ", round(pred, 2),
                                    "<br>Obs = ", round(obs_mean, 2), " +/- ", round(err, 2),
                                    "<br>Residual = ", round(residual, 2)),
                     hoverinfo = "text")
    }
    p %>% layout(xaxis  = list(title = paste0("Model predicted ", label), range = ax_range, gridcolor = "#eee"),
                yaxis  = list(title = paste0("Observed mean ", label),   range = ax_range, gridcolor = "#eee"),
                legend = list(orientation = "h", y = -0.28), margin = list(t = 5),
                plot_bgcolor = "#fefefe", paper_bgcolor = "#fefefe")
  })

  output[[paste0(prefix, "_plt_residuals")]] <- renderPlotly({
    res <- sim_results_fn(); req(!is.null(res))
    cd <- tryCatch(comp_data(), error = function(e) NULL)
    req(!is.null(cd) && nrow(cd) > 0)
    groups <- sort(unique(cd$group))
    p <- plot_ly()
    for (grp in groups) {
      d_g <- cd[cd$group == grp & !is.na(cd$residual), ]
      d_g <- d_g[order(d_g$time_h), ]
      col_vec <- ifelse(d_g$residual >= 0, "#2980b9", "#c0392b")
      lbl <- if (length(groups) <= 1) "Residual" else paste0("Residual: ", grp)
      p <- add_bars(p, data = d_g, x = ~time_h, y = ~residual, name = lbl,
                    marker = list(color = col_vec, line = list(color = "white", width = 0.5)),
                    text = ~paste0("t = ", round(time_h, 1), " h<br>Residual = ", sprintf("%+.2f", residual)),
                    hoverinfo = "text")
    }
    p %>% layout(xaxis = fr_xaxis_clock(res,fr_clock_step_from_input(input)),
                yaxis = list(title = "Residual (pred - obs)", zeroline = TRUE,
                             zerolinecolor = "#333", zerolinewidth = 1.5, gridcolor = "#eee"),
                barmode = "group", legend = list(orientation = "h", y = -0.35), margin = list(t = 5),
                plot_bgcolor = "#fefefe", paper_bgcolor = "#fefefe")
  })

  output[[paste0(prefix, "_tbl")]] <- renderDT({
    cd <- tryCatch(comp_data(), error = function(e) NULL)
    req(!is.null(cd) && nrow(cd) > 0)
    err_col <- if (isTRUE(input[[err_type_input]] == "sem")) "obs_sem" else "obs_sd"
    err_lbl <- if (isTRUE(input[[err_type_input]] == "sem")) "Obs SEM" else "Obs SD"
    disp <- cd %>% dplyr::arrange(time_h, group) %>%
      dplyr::transmute(
        `Time (h)` = round(time_h, 2),
        Group      = group,
        `Obs mean` = round(obs_mean, 2),
        !!err_lbl  := round(.data[[err_col]], 2),
        N          = obs_n,
        Model      = round(pred, 2),
        Residual   = round(residual, 2),
        `Abs error`= round(abs_error, 2)
      )
    datatable(disp, options = list(pageLength = 20, dom = "tp", scrollX = TRUE),
              rownames = FALSE, class = "compact stripe")
  })

  output[[paste0(prefix, "_dl")]] <- downloadHandler(
    filename = function() paste0("farima_comparison_", prefix, "_", Sys.Date(), ".csv"),
    content  = function(file) {
      cd <- tryCatch(comp_data(), error = function(e) NULL)
      req(!is.null(cd))
      write.csv(cd, file, row.names = FALSE)
    }
  )

  invisible(comp_data)
}

# Per-metric "suggest calibration" panel: lets the user check which published
# parameters the optimizer may search over, runs fr_suggest_calibration() on
# demand (never on the normal debounced auto-recalc path), and offers an
# Apply button that writes the suggestion back onto the real custom-parameter
# sliders. arm is which model arm ("PVT"/"KSS") supplies the parameters;
# comparison_metric is which fr_compare() bucket to score against -- these
# differ for the VAS tab, which is scored against the KSS prediction columns
# (fr_compare has no "VAS" case and falls through to the KSS default) but
# should still let the user search the KSS arm's parameters.
make_calibration_suggester <- function(input, output, session, prefix, arm, comparison_metric,
                                        param_choices, sim_inputs_fn) {
  rv_sugg <- reactiveVal(NULL)
  running <- reactiveVal(FALSE)
  shared_params <- c("eta_w","beta","zeta","phi","lambda","mu_s","nu","gamma")

  output[[paste0(prefix,"_calib_params_ui")]] <- renderUI({
    version <- input$model_version; req(version)
    avail <- param_choices
    if (version == "2021") avail <- setdiff(avail, c("gamma","nu"))
    if (version != "2024fe") avail <- setdiff(avail, "zeta")
    labels <- setNames(avail, vapply(avail, function(nm) {
      desc <- unname(FR_PARAM_INFO[nm])
      if (is.na(desc)) nm else paste0(nm, " — ", desc)
    }, character(1)))
    checkboxGroupInput(paste0(prefix,"_calib_params"), NULL, choices = labels, inline = FALSE)
  })

  output[[paste0(prefix,"_calib_warn")]] <- renderText({
    chosen <- input[[paste0(prefix,"_calib_params")]]
    if (length(chosen) > 6)
      "More than 6 free parameters on this few observations risks an ill-posed / overfit search."
    else ""
  })

  output[[paste0(prefix,"_calib_ready")]] <- reactive({
    x <- tryCatch(sim_inputs_fn(), error = function(e) NULL)
    isTRUE(input$custom_params) && !is.null(x) && isTRUE(x$ok) &&
      length(input[[paste0(prefix,"_calib_params")]]) > 0
  })
  outputOptions(output, paste0(prefix,"_calib_ready"), suspendWhenHidden = FALSE)

  observeEvent(input[[paste0(prefix,"_calib_run")]], {
    if (isTRUE(running())) return()
    chosen <- input[[paste0(prefix,"_calib_params")]]
    req(length(chosen) > 0)
    x <- isolate(sim_inputs_fn()); req(isTRUE(x$ok))
    running(TRUE)
    updateActionButton(session, paste0(prefix,"_calib_run"), label = "Optimizing...")
    res <- tryCatch(
      fr_suggest_calibration(x$value, arm, comparison_metric, chosen,
                              scale = isolate(input$prediction_scale)),
      error = function(e) list(error = conditionMessage(e)))
    running(FALSE)
    updateActionButton(session, paste0(prefix,"_calib_run"), label = "Suggest calibration")
    rv_sugg(res)
  })

  output[[paste0(prefix,"_calib_table")]] <- renderTable({
    r <- rv_sugg(); req(!is.null(r), is.null(r$error))
    data.frame(Parameter = r$param_names,
      Current = signif(r$before,4), Suggested = signif(r$after,4),
      `At bound` = ifelse(r$at_bound,"yes",""),
      `Also affects` = ifelse(r$param_names %in% shared_params,
                               if (arm=="PVT") "KSS (shared slider)" else "PVT (shared slider)", ""),
      check.names = FALSE)
  })
  output[[paste0(prefix,"_calib_ccc")]] <- renderText({
    r <- rv_sugg(); req(!is.null(r), is.null(r$error))
    sprintf("CCC: %.3f -> %.3f", r$ccc_before, r$ccc_after)
  })
  output[[paste0(prefix,"_calib_error")]] <- renderText({
    r <- rv_sugg(); if (!is.null(r$error)) r$error else ""
  })

  observeEvent(input[[paste0(prefix,"_calib_apply")]], {
    r <- rv_sugg(); req(!is.null(r), is.null(r$error), isTRUE(input$custom_params))
    for (i in seq_along(r$param_names)) {
      nm <- r$param_names[i]
      target_metric <- if (nm %in% shared_params) "PVT" else arm
      updateSliderInput(session, paste0("par_", target_metric, "_", nm), value = unname(r$after[i]))
    }
  })
  invisible(NULL)
}



# Use header delimiters for text imports; ISO timestamp colons are not separators.
classic_read_table <- function(path,filename) {
  if(tolower(tools::file_ext(filename)) %in% c("xlsx","xls","sav","zsav"))return(fr_read_table(path,filename))
  header<-readLines(path,n=1,warn=FALSE,encoding="UTF-8")
  if(!length(header))stop("The file is empty.",call.=FALSE)
  sep<-c(comma=",",semicolon=";",tab="\t")
  counts<-vapply(sep,function(mark){
    pattern<-paste0(mark,'(?=(?:[^"]*"[^"]*")*[^"]*$)')
    hits<-gregexpr(pattern,header,perl=TRUE)[[1]]
    sum(hits>0)
  },integer(1))
  chosen<-if(max(counts)>0)names(which.max(counts)) else "comma"
  fr_read_table(path,filename,chosen)
}

# Keep the preliminary app's automatic single-line time mapping.
# number_mode="index" (default): a bare number after the prefix (e.g. col_10)
# is the 10th occasion, spaced by fallback_step. number_mode="clock": that same
# bare number is a clock hour (e.g. col_10 = 10:00), needed for datasets like
# "Lapses_10, Lapses_12, ..., Lapses_2, Lapses_4" where the suffix wraps past
# midnight and isn't one of the hardcoded kss/vas/pvt column-name prefixes.
propose_times_from_names <- function(col_names,t0_clock=8,fallback_step=8,number_mode="index") {
  n<-length(col_names)
  if(!n)return(list(times=numeric(),method="none",detail="No columns selected."))
  # Loose gate: any digit anywhere in the name is enough to attempt the
  # clock-hour reading; fr_bulk_times()/fr_time_token() do the real,
  # more careful extraction and this call is tryCatch-guarded below, so a
  # permissive gate just means "try it," not "assume it will succeed."
  clock_form<-grepl("[0-9]",col_names)
  if(identical(number_mode,"clock") && all(clock_form)) {
    tt<-tryCatch(fr_bulk_times(col_names,"clock",t0_clock=t0_clock),error=function(e)NULL)
    if(!is.null(tt))return(list(times=unname(tt),method="clock-hour assumption",
      detail="Numbers after the underscore are treated as clock hours, in column order, across midnight. Text suffixes do not define time. Check this against the protocol."))
  }
  tt<-fr_propose_times(col_names,t0_clock,fallback_step)
  if(all(is.finite(tt)))return(list(times=unname(tt),method="names",
    detail="Times inferred from clock/day labels or occasion indices. Check the proposed times against the protocol."))
  list(times=seq(0,by=fallback_step,length.out=n),method="assumed regular spacing",
    detail=sprintf("Names are ambiguous: assuming first measurement at study start and %g-hour spacing in column order. This is a convenience assumption, not recovered timing; edit this one line if needed.",fallback_step))
}

classic_times <- function(text,n) {
  s<-trimws(text)
  if(!nzchar(s))stop("Enter measurement times in the time-mapping text box.",call.=FALSE)
  tokens<-strsplit(s,"[,;[:space:]]+")[[1]]
  v<-fr_parse_numbers(tokens,".","Measurement times",FALSE)
  if(length(v)!=n)stop("The number of measurement times must match the selected columns.",call.=FALSE)
  if(anyDuplicated(v))stop("Measurement times must be distinct; aggregate same-time replicates first.",call.=FALSE)
  v
}

make_wide_outcome_importer <- function(input,output,session,rec,rv,prefix,label,col_regex,default_step_h=2,origin_fn) {
  force(prefix);force(label);force(default_step_h)
  id<-function(s)paste0(prefix,"_",s)
  state<-reactiveValues(data=NULL,error=NULL,file=NULL,auto=TRUE,text="",cols=character(),proposal="",last_clock=NULL)
  observeEvent(input[[id("file")]],{
    f<-input[[id("file")]];req(f)
    state$data<-NULL;state$error<-NULL;state$file<-f;state$auto<-TRUE;state$cols<-character();state$text<-""
    tryCatch({state$data<-classic_read_table(f$datapath,f$name)},error=function(e){state$error<-conditionMessage(e)})
  },ignoreNULL=TRUE)
  observeEvent(input[[id("remove")]],{state$data<-NULL;state$error<-NULL;state$file<-NULL;state$text<-"";state$cols<-character()})
  output[[id("loaded")]]<-reactive(!is.null(state$data))
  outputOptions(output,id("loaded"),suspendWhenHidden=FALSE)
  output[[id("import_status")]]<-renderText({
    if(!is.null(state$error))return(paste("Import failed:",state$error))
    if(is.null(state$data))return("")
    map<-attr(state$data,"column_name_map");map<-map[map$changed,,drop=FALSE]
    if(!nrow(map))return("")
    paste("Renamed headings (all columns retained):",paste(paste0("column ",map$position," → ",map$imported),collapse="; "))
  })
  output[[id("ui_cols")]]<-renderUI({
    req(!is.null(state$data));cols<-names(state$data)
    auto<-grep(paste0("^",prefix,"([_0-9]|$)"),cols,value=TRUE,ignore.case=TRUE)
    pickerInput(id("cols_selected"),paste(label,"time-point columns:"),choices=cols,selected=auto,multiple=TRUE,
      options=list(`actions-box`=TRUE,`live-search`=TRUE,`selected-text-format`="count > 5",`count-selected-text`="{0} columns selected"))
  })
  output[[id("ui_id")]]<-renderUI({req(!is.null(state$data));selectInput(id("id_col"),"Subject ID column:",
    choices=c("(row number)"="__rownum__",names(state$data)),selected=names(state$data)[1])})
  output[[id("ui_group")]]<-renderUI({req(!is.null(state$data));selectInput(id("group_col"),"Group variable (optional):",
    choices=c("None"="__none__",names(state$data)),selected="__none__")})
  cols<-reactive({if(is.null(state$data))character() else intersect(input[[id("cols_selected")]],names(state$data))})
  # Auto-detect the raw scale's range from the selected columns' own data (any
  # number of categories: 0-500, 1-7, etc.), so vas_min/vas_max don't have to
  # be guessed or left at the 0-100 default. Mirrors the auto/manual toggle
  # already used for measurement-time detection below: a fresh detection
  # overwrites the fields only while state$range_auto is TRUE, and a real user
  # edit (value differing from what we last set) turns auto-detection off
  # until "Re-detect range" is clicked.
  if (prefix=="vas") {
    state$range_auto<-TRUE;state$range_lo<-NULL;state$range_hi<-NULL
    observe({
      req(!is.null(state$data));cc<-cols();req(length(cc)>0)
      dec<-if(is.null(input[[id("decimal")]]))"."else input[[id("decimal")]]
      v<-tryCatch(unlist(lapply(cc,function(cn)suppressWarnings(fr_parse_numbers(state$data[[cn]],dec,cn)))),
        error=function(e)numeric())
      v<-v[is.finite(v)];req(length(v)>0)
      lo<-floor(min(v));hi<-ceiling(max(v));if(hi<=lo)hi<-lo+1
      isolate({
        if(isTRUE(state$range_auto)){
          state$range_lo<-lo;state$range_hi<-hi
          updateNumericInput(session,"vas_min",value=lo)
          updateNumericInput(session,"vas_max",value=hi)
        }
      })
    },priority=9)
    observeEvent(input$vas_min,{if(!identical(input$vas_min,state$range_lo))state$range_auto<-FALSE},ignoreInit=TRUE)
    observeEvent(input$vas_max,{if(!identical(input$vas_max,state$range_hi))state$range_auto<-FALSE},ignoreInit=TRUE)
    observeEvent(input$vas_redetect_range,{state$range_auto<-TRUE})
    output$vas_range_note<-renderText({
      req(!is.null(state$data));cc<-cols();req(length(cc)>0)
      if(isTRUE(state$range_auto))sprintf("Auto-detected from %d column(s): %s to %s.",length(cc),
        format(state$range_lo),format(state$range_hi))
      else "Manually edited. Click \"Re-detect range\" to replace with the observed data range."
    })
  }
  proposal<-reactive({propose_times_from_names(cols(),origin_fn()$t0_clock,
    if(is.null(input[[id("step")]]))default_step_h else input[[id("step")]],
    if(is.null(input[[id("time_number_mode")]]))"index" else input[[id("time_number_mode")]])})
  # The original compact input remains; no table and no per-cell entry.
  output[[id("ui_times")]]<-renderUI({
    req(!is.null(state$data))
    isolate(tagList(
      textInput(id("times_txt"),"Measurement times (h from t=0), same order as columns:",state$text,placeholder="e.g. 0, 8, 16, 24, 32"),
      textOutput(id("time_note")),
      radioButtons(id("time_number_mode"),"A bare number after the column prefix (e.g. name_10) means:",
        choices=c("Sequential occasion index"="index","Clock hour of day"="clock"),selected="index",inline=TRUE),
      numericInput(id("step"),"Spacing if names have no clock time (h):",default_step_h,min=.001),
      actionButton(id("btn_redetect"),"Re-detect from column names",class="btn-xs btn-default",icon=icon("sync"))))
  })
  observe({
    req(!is.null(state$data));cc<-cols();prop<-proposal();clock<-origin_fn()$t0_clock
    isolate({
      changed<-!identical(cc,state$cols)
      if(changed||state$auto){
        state$auto<-TRUE;state$text<-paste(prop$times,collapse=", ");state$proposal<-prop$detail
        updateTextInput(session,id("times_txt"),value=state$text)
      }
      state$cols<-cc;state$last_clock<-clock
    })
  },priority=10)
  observeEvent(input[[id("times_txt")]],{
    txt<-input[[id("times_txt")]]
    if(!identical(txt,state$text)) {state$text<-txt;state$auto<-FALSE}
  },ignoreInit=TRUE)
  observeEvent(input[[id("btn_redetect")]],{
    p<-proposal();state$text<-paste(p$times,collapse=", ");state$proposal<-p$detail;state$auto<-TRUE
    updateTextInput(session,id("times_txt"),value=state$text)
  })
  output[[id("time_note")]]<-renderText({req(!is.null(state$data));
    paste(length(cols()),"columns.",if(state$auto)paste("Auto-detected:",state$proposal) else "Manual time mapping retained. Re-detect to replace it from column names.")})
  result<-reactive({
    if(!is.null(state$error))stop(label," import: ",state$error,call.=FALSE)
    d<-state$data;if(is.null(d))return(NULL)
    cc<-cols();if(!length(cc))return(list(observations=NULL,settings=list(columns=character()),
      source=list(filename=state$file$name,raw_table=d),disabled=TRUE))
    times<-classic_times(state$text,length(cc))
    subj<-input[[id("id_col")]];grp<-input[[id("group_col")]]
    req(!is.null(subj),!is.null(grp))
    settings<-list(columns=cc,id=if(subj=="__rownum__")"__row__" else subj,group=grp,times=times,
      decimal=input[[id("decimal")]],vas_min=input$vas_min,vas_max=input$vas_max,direction=input$vas_direction,
      target=if(is.null(input$vas_target))"kss"else input$vas_target,
      pvt_units_confirmed=TRUE,lapse_rule=if(is.null(input$pvt_lapse_rule))">=500 ms"else input$pvt_lapse_rule)
    obs<-fr_prepare_observations(d,settings,label)
    settings$time_mapping<-if(state$auto)state$proposal else "User-entered elapsed hours"
    list(observations=obs,settings=settings,disabled=FALSE,
      source=list(filename=state$file$name,sha256=digest::digest(file=state$file$datapath,algo="sha256"),
        raw_table=d,column_name_map=attr(d,"column_name_map")))
  })
  output[[id("summary_txt")]]<-renderText({
    tryCatch({r<-result();if(is.null(r))return("");if(isTRUE(r$disabled))return("No outcome columns selected; this dataset is excluded.")
      d<-r$observations;v<-d$value_raw;v<-v[is.finite(v)]
      sprintf("Subjects: %d | Columns: %d\nValid observations: %d | Missing: %d\nMean: %s | SD: %s",
        length(unique(d$subject)),length(r$settings$columns),length(v),sum(is.na(d$value_raw)),
        if(length(v))format(mean(v),digits=4)else"--",if(length(v)>1)format(sd(v),digits=4)else"--")
    },error=function(e)conditionMessage(e))
  })
  output[[id("tbl_preview")]]<-renderDT({req(!is.null(state$data));
    datatable(head(state$data,10),options=list(pageLength=5,dom="t",scrollX=TRUE),rownames=FALSE,class="compact stripe")})
  result
}

# Row of fit-statistic value boxes for one comparison sub-tab (KSS/VAS/PVT/
# Index), each with a short explanation of what it measures underneath so the
# meaning doesn't have to be looked up elsewhere. Server side is wired by
# make_comparison_outputs()'s _vb_rmse/_vb_mae/_vb_r/_vb_ccc/_vb_ccc_ind
# outputs.
# A collapsed <details>/<summary> disclosure per stat -- no server-side
# wiring needed, hidden until clicked. `question` is the short framing line
# shown always; `body` (a tagList of paragraphs) only renders when expanded.
fr_stat_explainer <- function(question, ...) {
  tags$details(style = "margin-top: 2px; font-size: 85%; color: #555;",
    tags$summary(style = "cursor: pointer; color: #337ab7;", question),
    tags$div(style = "margin-top: 4px;", ...)
  )
}

fr_fit_stat_boxes <- function(prefix) {
  fluidRow(
    column(2, valueBoxOutput(paste0(prefix,"_vb_rmse"), width = 12),
      fr_stat_explainer("Prediction error, with larger errors weighted more heavily",
        "Every miss is squared before averaging, then the result is square-rooted back into the outcome's own units.",
        tags$br(), "Example: if the model misses by 1 point on one occasion and by 5 points on another, RMSE weighs that 5-point miss much more heavily than the 1-point one -- so it's pulled up by a few bad misses more than by many small ones.")),
    column(2, valueBoxOutput(paste0(prefix,"_vb_mae"), width = 12),
      fr_stat_explainer("Typical size of the prediction error",
        "Same idea as RMSE, but every miss counts in direct proportion to its size -- no squaring.",
        tags$br(), "Example: misses of 1 and 5 points average to an MAE of 3, right in the middle. It's usually a bit smaller than RMSE, and doesn't get dragged up by the occasional large miss.")),
    column(2, valueBoxOutput(paste0(prefix,"_vb_r"), width = 12),
      fr_stat_explainer("Similarity of the predicted and observed pattern",
        "This only checks whether predicted and observed values rise and fall together -- not whether they're actually close.",
        tags$br(), "Example: if the observed curve goes 3 -> 5 -> 7 and the model predicts 6 -> 10 -> 14, the two move in perfect lockstep, so r can sit right near 1 -- even though every prediction is about double the true value.",
        tags$br(), "r measures similarity of pattern, not agreement in size.")),
    column(2, valueBoxOutput(paste0(prefix,"_vb_ccc"), width = 12),
      fr_stat_explainer("Agreement of the predicted and observed average curves",
        "First, every participant's value at each measurement occasion is averaged together. Then the predicted curve is compared against that observed ",tags$i("average"),
        " curve -- checking both whether they move together (like r) and whether they actually land close together (unlike r).",
        tags$br(), "A high value means the model reproduces the average trajectory well -- but averaging across people first also removes most of the person-to-person noise, so this can look good even if it wouldn't for any one individual.")),
    column(2, valueBoxOutput(paste0(prefix,"_vb_ccc_ind"), width = 12),
      fr_stat_explainer("Agreement of predictions with individual measurements",
        "Now nobody is averaged away: each person's own raw measurement is compared directly with what the model predicted for that occasion -- a much harder test than the group version.",
        tags$br(), "A big drop from CCC (group) to here usually means the model captures the population-average pattern well, but individual people deviate substantially from that average.")),
    column(2, valueBoxOutput(paste0(prefix,"_vb_ccc_subj"), width = 12),
      fr_stat_explainer("Average agreement within each person over time",
        "Here CCC is computed separately within each person -- using only that person's own repeated measurements -- and those personal scores are then averaged across people.",
        tags$br(), "This asks something different from CCC (individual) above: not \"how close is each point to the model\" but \"does the shape of one person's own curve follow the model's shape over time,\" independent of that person's overall level.",
        tags$br(), "Subjects with fewer than 2 usable points, or a perfectly flat response, can't have a meaningful personal score and are dropped (see n used/total).",
        tags$br(), tags$b("Caveat: "),"this is a simple average of independent per-subject estimates. Tsai et al. (2017, Statistics in Medicine, doi:10.1002/sim.7207) describe GEE- and variance-components-based estimators that instead model the correlation between a subject's repeated measurements directly, which is the more statistically rigorous way to pool this kind of agreement -- not implemented here."))
  )
}

# UI block for the "Suggest calibration" panel on the KSS/VAS/PVT comparison
# sub-tabs. Hidden unless "Explore custom parameters" is on, since Apply
# writes into the custom-parameter sliders, which don't exist otherwise.
fr_calib_box <- function(prefix) {
  tagList(
    conditionalPanel("input.custom_params",
      fluidRow(box(collapsible=TRUE,title = "Suggest calibration", width = 12, status = "info", solidHeader = TRUE,
        helpText("Searches only the checked parameters, within the same bounds as their sliders, to maximize CCC (Lin's concordance correlation coefficient) against the loaded, currently selected observations."),
        uiOutput(paste0(prefix, "_calib_params_ui")),
        textOutput(paste0(prefix, "_calib_warn")),
        actionButton(paste0(prefix, "_calib_run"), "Suggest calibration", class = "btn-info"),
        br(), br(),
        textOutput(paste0(prefix, "_calib_error")),
        tableOutput(paste0(prefix, "_calib_table")),
        textOutput(paste0(prefix, "_calib_ccc")),
        actionButton(paste0(prefix, "_calib_apply"), "Apply suggested values", class = "btn-success")
      ))
    ),
    conditionalPanel("!input.custom_params",
      fluidRow(box(collapsible=TRUE,width = 12, status = "warning",
        helpText('Enable "Explore custom parameters" on the Parameters tab to get calibration suggestions.'))))
  )
}

ui <- dashboardPage(
  skin = "blue",
  dashboardHeader(title = "FaRiMa 1.2"),
  dashboardSidebar(
    sidebarMenu(id = "tabs",
      menuItem("Schedule",   tabName = "schedule", icon = icon("bed")),
      menuItem("Parameters", tabName = "params",   icon = icon("sliders-h")),
      menuItem("Results",    tabName = "results",  icon = icon("chart-line")),
      menuItem("Comparison", tabName = "compare",  icon = icon("not-equal")),
      menuItem("Export",     tabName = "export",   icon = icon("download"))
    ),
    hr(),
    div(style = "padding: 10px;",
      p(style = "color:#aaa; font-size:11px;",
        "FaRiMa 1.2 -- corrected original interface",
        br(), "Published McCauley parameter sets (2021 / 2024)",
        br(), "Live recalculation: results update as you type/slide.")
    )
  ),

  dashboardBody(
    tags$head(
      tags$style(HTML("
        .content-wrapper { background-color: #f4f6f9; }
        .box-header .box-title { font-weight: 600; }
        pre { font-size: 11px; }
        .plotly.html-widget { resize: vertical; overflow: hidden; min-height: 150px; }
      ")),
      tags$script(HTML("
        (function(){
          function attach(el){
            if (el._farimaResizeObserved) return;
            el._farimaResizeObserved = true;
            var ro = new ResizeObserver(function(){
              if (window.Plotly) Plotly.Plots.resize(el);
            });
            ro.observe(el);
          }
          $(document).on('shiny:value shiny:visualchange', function(){
            document.querySelectorAll('.plotly.html-widget').forEach(attach);
          });
        })();
      "))
    ),

    tabItems(

      # ========== SCHEDULE TAB ==========
      tabItem(tabName = "schedule",
        fluidRow(
          box(collapsible=TRUE,title = "Study reference time", width = 4, status = "primary", solidHeader = TRUE,
            numericInput("t0_clock", "Clock hour at t = 0 (from midnight):",
                        value = 8, min = 0, max = 23.5, step = 0.5),
            numericInput("duration","Study duration (hours):",24,min=.01,step=1),
            checkboxInput("dated","Use dated schedule / actigraphy",FALSE),
            conditionalPanel("input.dated",textInput("start_iso","Study start (ISO, including UTC offset):","",placeholder="2026-09-18T08:00:00+02:00")),
            helpText("Sleep and observations share this origin. An empty schedule means continuous wake. For dated input, the starting clock comes from the entered datetime.",
              br(),"Auto-extends to cover the latest loaded measurement time + 0.5 h once KSS/VAS/PVT data is loaded; edit it manually to override.")
          ),
          box(collapsible=TRUE,title = "Input mode", width = 4, status = "primary", solidHeader = TRUE,
            radioButtons("input_mode", NULL,
              choices = c("Manual schedule entry" = "manual",
                          "Upload file (onset / offset)"  = "file",
                          "Upload actigraphy (datetime + state)" = "actig"),
              selected = "manual"),
            conditionalPanel("input.input_mode == 'file'",
              fileInput("schedule_file", "File (CSV / xlsx):",
                       accept = c(".csv", ".txt", ".xlsx", ".xls")),
              selectInput("schedule_format","Time format:",c("Elapsed hours"="hours","ISO datetime with UTC offset"="iso")),
              helpText("First two columns: onset and offset. Sleep opportunities are the published model input.")
            ),
            conditionalPanel("input.input_mode == 'actig'",
              fileInput("actig_file", "Actigraphy file (CSV):",
                       accept = c(".csv", ".txt")),
              numericInput("epoch_minutes","Epoch duration (minutes):",1,min=.001),
              helpText("Experimental: actual-sleep switching is not equivalent to the published sleep-opportunity calibration. Columns: datetime (ISO + offset), sleep_wake (0 sleep / 1 wake). Epochs must cover the full study.")
            )
          ),
          box(collapsible=TRUE,title = "Run", width = 4, status = "success", solidHeader = TRUE,
            actionButton("btn_example", "Load 7-day example",
                        class = "btn-default btn-block", icon = icon("magic")),
            br(),
            checkboxInput("auto_recalc",
                         "Live recalculation (auto-update on any change)",
                         value = TRUE),
            conditionalPanel("!input.auto_recalc",
              actionButton("btn_run", "Run Simulation",
                          class = "btn-success btn-lg btn-block")
            ),
            br(),
            verbatimTextOutput("sim_status_txt", placeholder = TRUE)
          )
        ),

        # Manual entry panel
        conditionalPanel("input.input_mode == 'manual'",
          fluidRow(
            box(collapsible=TRUE,title = "Sleep schedule (hours from t = 0)", width = 12,
                status = "info", solidHeader = TRUE,
              fluidRow(
                column(3, numericInput("new_onset", "Sleep onset (h from t=0):",
                                      value = 15, min = 0, max = 500, step = 0.5)),
                column(3, numericInput("new_offset", "Wake time (h from t=0):",
                                      value = 23, min = 0, max = 500, step = 0.5)),
                column(3, br(),
                  actionButton("btn_add", "Add period", class = "btn-info"),
                  actionButton("btn_del", "Delete selected", class = "btn-warning")
                ),
                column(3, br(),
                  helpText("Onset 15 h -> 23:00 if study starts at 08:00.",
                          br(), "Select row(s) to delete.")
                )
              ),
              DTOutput("tbl_schedule_manual")
            )
          )
        ),

        # File preview
        conditionalPanel("input.input_mode != 'manual'",
          fluidRow(
            box(collapsible=TRUE,title = "Parsed schedule preview (hours from t=0)", width = 12,
                status = "info", solidHeader = TRUE,
              DTOutput("tbl_file_preview")
            )
          )
        ),

        # ---- Empirical KSS data import (optional) -------------------
        fluidRow(
          box(collapsible=TRUE,title = "Empirical KSS overlay (optional)",
              width = 12, status = "warning", solidHeader = TRUE,
              collapsed = FALSE,
            fluidRow(
              column(4,
                fileInput("kss_file", "Upload KSS file (CSV / xlsx):",
                         accept = c(".csv", ".txt", ".xlsx", ".xls")),
                radioButtons("kss_decimal","Decimal mark:",c("Point"=".","Comma"=","),inline=TRUE),
                actionButton("kss_remove","Remove dataset",class="btn-xs btn-default"),
                verbatimTextOutput("kss_import_status"),
                helpText(
                  "Wide format: one row per subject, one column per KSS measurement.",
                  br(), "Column names matching kss1, kss2, ... are auto-detected.",
                  br(), "UTF-8 with BOM, comma- and semicolon-delimited files are supported."
                )
              ),
              column(4,
                conditionalPanel("output.kss_loaded",
                  uiOutput("kss_ui_cols"),
                  uiOutput("kss_ui_id"),
                  uiOutput("kss_ui_group")
                )
              ),
              column(4,
                conditionalPanel("output.kss_loaded",
                  uiOutput("kss_ui_times"),
                  br(),
                  checkboxInput("kss_show_indiv", "Show individual subject traces", value = FALSE),
                  checkboxInput("kss_show_mean",  "Show mean +/- SD band", value = TRUE)
                )
              )
            ),
            conditionalPanel("output.kss_loaded",
              hr(),
              fluidRow(
                column(5, verbatimTextOutput("kss_summary_txt")),
                column(7, DTOutput("kss_tbl_preview"))
              )
            )
          )
        ),

        # ---- Empirical VAS data import (optional) --------------------
        fluidRow(
          box(collapsible=TRUE,title = "Empirical VAS / other-scale overlay (optional, any numeric scale)",
              width = 12, status = "warning", solidHeader = TRUE,
              collapsed = TRUE,
            fluidRow(
              column(4,
                fileInput("vas_file", "Upload VAS / other-scale file (CSV / xlsx):",
                         accept = c(".csv", ".txt", ".xlsx", ".xls")),
                radioButtons("vas_decimal","Decimal mark:",c("Point"=".","Comma"=","),inline=TRUE),
                actionButton("vas_remove","Remove dataset",class="btn-xs btn-default"),
                verbatimTextOutput("vas_import_status"),
                helpText(
                  "Wide format, one column per measurement. Works for any single ordinal/continuous",
                  " sleepiness-type scale, not just a 0-100 VAS (e.g. a 1-7 or 0-500 instrument).",
                  br(), "Columns matching vas1, vas2, ... are auto-detected;",
                  " otherwise select manually below (e.g. custom-named columns)."
                )
              ),
              column(4,
                conditionalPanel("output.vas_loaded",
                  uiOutput("vas_ui_cols"),
                  uiOutput("vas_ui_id"),
                  uiOutput("vas_ui_group")
                )
              ),
              column(4,
                conditionalPanel("output.vas_loaded",
                  uiOutput("vas_ui_times"),
                  br(),
                  checkboxInput("vas_show_indiv", "Show individual subject traces", value = FALSE),
                  checkboxInput("vas_show_mean",  "Show mean +/- SD band", value = TRUE)
                )
              )
            ),
            conditionalPanel("output.vas_loaded",
              hr(),
              fluidRow(
                column(3, numericInput("vas_min", "Raw scale minimum:", value = 0,   step = 1)),
                column(3, numericInput("vas_max", "Raw scale maximum:", value = 100, step = 1)),
                column(3, radioButtons("vas_direction", "Scale direction:",
                                       choices = c("Low end = alert, high end = sleepy" = "direct",
                                                   "Low end = sleepy, high end = alert" = "reverse"),
                                       selected = "direct")),
                column(3, radioButtons("vas_target", "Rescale to:",
                                       choices = c("KSS-equivalent (1–9)" = "kss",
                                                   "0–1 (compare to composite risk index)" = "01"),
                                       selected = "kss"))
              ),
              fluidRow(
                column(9, textOutput("vas_range_note")),
                column(3, actionButton("vas_redetect_range","Re-detect range",class="btn-xs btn-default",icon=icon("sync")))
              ),
              helpText(
                "NOT VALIDATED: the raw scale is rescaled linearly to the target above for exploratory overlays; this is not a validated conversion. Min/max are auto-detected from the observed data and can be edited",
                " if the instrument's true endpoints differ from what happened to be observed in this sample. Values outside the stated endpoints are rejected, not silently clipped."
              ),
              fluidRow(
                column(5, verbatimTextOutput("vas_summary_txt")),
                column(7, DTOutput("vas_tbl_preview"))
              )
            )
          )
        ),

        # ---- Empirical PVT data import (optional) --------------------
        fluidRow(
          box(collapsible=TRUE,title = "Empirical PVT overlay (optional)",
              width = 12, status = "warning", solidHeader = TRUE,
              collapsed = TRUE,
            fluidRow(
              column(4,
                fileInput("pvt_file", "Upload PVT lapse file (CSV / xlsx):",
                         accept = c(".csv", ".txt", ".xlsx", ".xls")),
                radioButtons("pvt_decimal","Decimal mark:",c("Point"=".","Comma"=","),inline=TRUE),
                actionButton("pvt_remove","Remove dataset",class="btn-xs btn-default"),
                selectInput("pvt_lapse_rule","Recorded PVT lapse definition:",c(">=500 ms",">500 ms")),
                verbatimTextOutput("pvt_import_status"),
                helpText(
                  "Expected metric: ", tags$b("lapses (count per 10-min bout)"),
                  " -- the same unit as the model's PVT arm.",
                  br(), "Columns matching pvt1, pvt2, ... are auto-detected.",
                  br(), tags$em("Reaction time is not accepted here: the model predicts lapse",
                                " counts, not latencies, so RT has no comparable model arm.")
                )
              ),
              column(4,
                conditionalPanel("output.pvt_loaded",
                  uiOutput("pvt_ui_cols"),
                  uiOutput("pvt_ui_id"),
                  uiOutput("pvt_ui_group")
                )
              ),
              column(4,
                conditionalPanel("output.pvt_loaded",
                  uiOutput("pvt_ui_times"),
                  br(),
                  checkboxInput("pvt_show_indiv", "Show individual subject traces", value = FALSE),
                  checkboxInput("pvt_show_mean",  "Show mean +/- SD band", value = TRUE)
                )
              )
            ),
            conditionalPanel("output.pvt_loaded",
              hr(),
              fluidRow(
                column(5, verbatimTextOutput("pvt_summary_txt")),
                column(7, DTOutput("pvt_tbl_preview"))
              )
            )
          )
        )

      ),  # end schedule tabItem


      # ========== PARAMETERS TAB ==========
      tabItem(tabName = "params",
        fluidRow(
          box(collapsible=TRUE,title="Published model",width=3,status="primary",solidHeader=TRUE,
            selectInput("model_version","Parameter set:",c("2024 recovery + inertia"="2024fe","2024 sleep inertia"="2024jtb","2021 sleepiness"="2021")),
            helpText("Defaults reproduce published equations and calibrated group-level parameters. This app has numerical verification; it is not independently validated for your population."),
            uiOutput("model_reference")),
          box(collapsible=TRUE,title="Initial conditions",width=3,status="primary",solidHeader=TRUE,
            numericInput("baseline_on","Habitual sleep onset (clock h):",23,min=0,max=23.99,step=.5),
            numericInput("baseline_off","Habitual awakening (clock h):",7,min=0,max=23.99,step=.5),
            helpText("A stable repeating baseline is aligned to study start. This is a pre-study sleep-history assumption, not a measured individual state.")),
          box(collapsible=TRUE,title="Simulation settings",width=3,status="warning",solidHeader=TRUE,
            selectInput("dt_choice","Output spacing:",c("1 minute (default)"="1","5 minutes"="5","15 minutes"="15"),selected="1"),
            radioButtons("prediction_scale","Display / comparison:",c("Raw predictions"="raw","Bounded to outcome scale"="bounded"),selected="raw"),
            helpText("Adaptive integration uses exact sleep/wake boundaries. Output spacing affects plots and exposure precision.")),
          box(collapsible=TRUE,title="Calibration status",width=3,status="warning",solidHeader=TRUE,
            checkboxInput("custom_params","Explore custom parameters",FALSE),
            helpText("Published values are used unless this is enabled. Edited values and phase shifts are unvalidated sensitivity analyses. No automatic chronotype-to-phase conversion is claimed."),
            textOutput("calibration_status"))
        ),
        fluidRow(
          box(collapsible=TRUE,title="Homeostatic processes — published rates",width=4,status="warning",solidHeader=TRUE,
            uiOutput("homeostatic_controls")),
          box(collapsible=TRUE,title="Allostatic process — published feedback",width=4,status="danger",solidHeader=TRUE,
            uiOutput("allostatic_controls")),
          box(collapsible=TRUE,title="Circadian process — published sine model",width=4,status="info",solidHeader=TRUE,
            uiOutput("circadian_controls"))
        ),
        fluidRow(
          box(collapsible=TRUE,title="Sleep inertia — published rate and scaling",width=4,status="primary",solidHeader=TRUE,
            uiOutput("inertia_controls")),
          box(collapsible=TRUE,title="KSS output — published scaling",width=4,status="success",solidHeader=TRUE,
            uiOutput("kss_parameter_controls")),
          box(collapsible=TRUE,title="PVT output — published scaling",width=4,status="success",solidHeader=TRUE,
            uiOutput("pvt_parameter_controls"))
        ),
        fluidRow(
          box(collapsible=TRUE,title="Composite impairment index — NOT VALIDATED",width=6,status="danger",solidHeader=TRUE,
            sliderInput("w1","KSS weight:",0,1,.35,.05),
            sliderInput("w2","PVT weight:",0,1,.40,.05),
            numericInput("pvt_ceiling","PVT normalization ceiling (lapses / 10 min):",20,min=.01),
            numericInput("index_threshold","Illustrative index threshold:",.66,min=0,max=1,step=.01),
            textOutput("weight_check"),
            helpText("Author-defined index, not accident probability. Weights, the ceiling of 20 and threshold 0.66 have no demonstrated validation. Weights are normalized to sum to one. The preliminary S_eff term was removed because it belonged to the unsupported model."),
            tags$code("Index = a × (bounded KSS − 1)/8 + b × min(max(PVT,0)/ceiling,1)")),
          box(collapsible=TRUE,title="What has evidence, and what is a choice?",width=6,status="info",solidHeader=TRUE,
            tags$ul(tags$li("Published calibration: separate PVT/KSS equations, constants and sleep-inertia rate."),
              tags$li("Study assumptions: schedule, duration, pre-study sleep history and measurement times."),
              tags$li("Unvalidated extensions: composite index, its weights/threshold, VAS-to-1–9 rescaling, custom model parameters and actual-sleep actigraphy switching.")),
            helpText("Numerical agreement with a publication does not establish predictive accuracy in a new person, population or task. Comparison panels describe fit to your uploaded data; they do not fit or validate coefficients."))
        )
      ),

      # ========== RESULTS TAB ==========
      tabItem(tabName = "results",
        conditionalPanel("!output.sim_done",
          fluidRow(box(collapsible=TRUE,width = 12, status = "warning",
            h4("Set a valid study duration and schedule on the Schedule tab; an empty schedule represents continuous wake.")
          ))
        ),
        conditionalPanel("output.sim_done",
          fluidRow(
            valueBoxOutput("vb_kss",      width = 4),
            conditionalPanel("output.pvt_loaded", valueBoxOutput("vb_pvt", width = 4)),
            valueBoxOutput("vb_highrisk", width = 4)
          ),
          fluidRow(
            box(collapsible=TRUE,width = 12, status = "primary", solidHeader = FALSE,
              selectInput("clock_tick_hours", "Time axis resolution:",
                choices = c("Auto"="", "Every 1 h"="1", "Every 2 h"="2", "Every 3 h"="3",
                            "Every 4 h"="4", "Every 6 h"="6", "Every 8 h"="8",
                            "Every 12 h"="12", "Every 24 h"="24"),
                selected = "", width = "220px"),
              helpText("Controls the tick spacing on every time-series plot's clock-time axis (Results and Comparison tabs). Drag a plot's bottom-right corner to resize it.")
            )
          ),
          fluidRow(
            box(collapsible=TRUE,title = "Subjective sleepiness (KSS) & PVT lapses / 10 min — with empirical overlays",
                width = 12, status = "primary", solidHeader = TRUE,
              plotlyOutput("plt_kss_pvt", height = "340px"),
              helpText("Circles = observed KSS, triangles = VAS/other-scale rescaled to a KSS-equivalent, squares = observed PVT lapses (when loaded).")
            )
          ),
          fluidRow(
            box(collapsible=TRUE,title = "Composite impairment index — author-defined, NOT VALIDATED",
                width = 12, status = "danger", solidHeader = TRUE,
              fluidRow(
                column(5,
                  checkboxInput("risk_overlay_sleep",
                                "Overlay observed sleepiness (KSS / VAS) on right axis",
                                value = FALSE)),
                column(7,
                  conditionalPanel("input.risk_overlay_sleep",
                    radioButtons("risk_overlay_scale", NULL, inline = TRUE,
                      choices = c("Right axis in KSS units (1-9)" = "kss",
                                  "Mapped onto the index scale [0-1]" = "risk"),
                      selected = "kss")))
              ),
              checkboxInput("risk_overlay_pvt",
                            "Overlay observed PVT (mapped onto the index scale [0-1], same left axis)",
                            value = FALSE),
              plotlyOutput("plt_risk", height = "280px"),
              helpText("Weights and threshold are author choices, not calibrated risk limits. The line marks the selected illustrative threshold; the index is not a probability."),
              conditionalPanel("input.risk_overlay_sleep && input.risk_overlay_scale == 'risk'",
                helpText("Mapped values use only the KSS term of the composite, (KSS-1)/8, so they are NOT an observed composite score -- the PVT term cannot be derived from a sleepiness rating."))
            )
          ),
          fluidRow(
            box(collapsible=TRUE,title = "State variables — published model",
                width = 12, status = "info", solidHeader = TRUE,
              selectInput("state_metric","Outcome:",c("PVT","KSS")),
              selectInput("state_variable","State:",c("Homeostatic q (circadian-coupled, not monotonic)"="q","Allostatic u"="u","Circadian amplitude kappa"="kappa","Circadian sine c"="c","Sleep inertia h"="h")),
              plotlyOutput("plt_states", height = "300px")
            )
          ),
          fluidRow(
            box(collapsible=TRUE,title = "Hourly prediction table (wake periods only)",
                width = 12, status = "success", solidHeader = TRUE,
              DTOutput("tbl_preds"),
              br(),
              downloadButton("dl_csv", "Download full predictions CSV",
                            class = "btn-sm btn-success")
            )
          )
        )
      ),  # end results tabItem


      # ========== COMPARISON TAB ==========
      tabItem(tabName = "compare",

        conditionalPanel("!output.sim_done || (!output.kss_loaded && !output.vas_loaded && !output.pvt_loaded)",
          fluidRow(
            box(collapsible=TRUE,width = 12, status = "warning", solidHeader = FALSE,
              h4("Load at least one empirical dataset to activate this tab:"),
              tags$ol(
                tags$li("A valid schedule/simulation (Schedule tab)"),
                tags$li("KSS, VAS, and/or PVT data (Schedule tab, collapsible boxes below the schedule)")
              )
            )
          )
        ),

        conditionalPanel("output.sim_done && (output.kss_loaded || output.vas_loaded || output.pvt_loaded)",
          fluidRow(
            box(collapsible=TRUE,width = 12, status = "primary", solidHeader = FALSE,
              radioButtons("err_type", "Error bars show:",
                           choices  = c("+/- SD" = "sd", "+/- SEM" = "sem"),
                           selected = "sd", inline = TRUE),
              helpText("SD = individual variability; SEM = precision of the mean. Applies to all sub-tabs below.")
            )
          ),

          uiOutput("cross_modality_section"),

          tabsetPanel(id = "compare_subtabs",

            # ---- KSS sub-tab ----
            tabPanel("KSS",
              conditionalPanel("!output.kss_loaded",
                fluidRow(box(collapsible=TRUE,width = 12, status = "warning", h5("No empirical KSS data loaded.")))
              ),
              conditionalPanel("output.kss_loaded",
                fluidRow(
                  box(collapsible=TRUE,width = 12, status = "primary", solidHeader = FALSE,
                      uiOutput("kss_ui_groups"),uiOutput("kss_ui_subjects"),verbatimTextOutput("kss_exclusions"))
                ),
                fr_fit_stat_boxes("kss"),
                fr_calib_box("kss"),
                fluidRow(
                  box(collapsible=TRUE,title = "Model vs Observed KSS -- time series", width = 8,
                      status = "primary", solidHeader = TRUE,
                    plotlyOutput("kss_plt_ts", height = "360px")),
                  box(collapsible=TRUE,title = "Predicted vs Observed", width = 4,
                      status = "warning", solidHeader = TRUE,
                    plotlyOutput("kss_plt_scatter", height = "360px"))
                ),
                fluidRow(
                  box(collapsible=TRUE,title = "Residuals (Model - Observed)", width = 12,
                      status = "danger", solidHeader = TRUE,
                    plotlyOutput("kss_plt_residuals", height = "200px"))
                ),
                fluidRow(
                  box(collapsible=TRUE,title = "Comparison table", width = 12,
                      status = "success", solidHeader = TRUE,
                    DTOutput("kss_tbl"), br(),
                    downloadButton("kss_dl", "Download CSV", class = "btn-sm btn-success"))
                )
              )
            ),

            # ---- VAS sub-tab ----
            tabPanel("VAS (-> KSS-rescaled)",
              conditionalPanel("!output.vas_loaded",
                fluidRow(box(collapsible=TRUE,width = 12, status = "warning", h5("No empirical VAS data loaded.")))
              ),
              conditionalPanel("output.vas_loaded",
                fluidRow(
                  box(collapsible=TRUE,width = 12, status = "primary", solidHeader = FALSE,
                      uiOutput("vas_ui_groups"),uiOutput("vas_ui_subjects"),verbatimTextOutput("vas_exclusions"))
                ),
                fr_fit_stat_boxes("vas"),
                fr_calib_box("vas"),
                fluidRow(
                  box(collapsible=TRUE,title = "Model vs VAS / other-scale overlay -- time series", width = 8,
                      status = "primary", solidHeader = TRUE,
                    plotlyOutput("vas_plt_ts", height = "360px")),
                  box(collapsible=TRUE,title = "Predicted vs Observed", width = 4,
                      status = "warning", solidHeader = TRUE,
                    plotlyOutput("vas_plt_scatter", height = "360px"))
                ),
                fluidRow(
                  box(collapsible=TRUE,title = "Residuals (Model - Observed)", width = 12,
                      status = "danger", solidHeader = TRUE,
                    plotlyOutput("vas_plt_residuals", height = "200px"))
                ),
                fluidRow(
                  box(collapsible=TRUE,title = "Comparison table", width = 12,
                      status = "success", solidHeader = TRUE,
                    DTOutput("vas_tbl"), br(),
                    downloadButton("vas_dl", "Download CSV", class = "btn-sm btn-success"))
                ),
                helpText("Rescaled via linear endpoint mapping of the raw scale to the chosen target (see Schedule tab note); treat as approximate.")
              )
            ),

            # ---- PVT sub-tab ----
            tabPanel("PVT",
              conditionalPanel("!output.pvt_loaded",
                fluidRow(box(collapsible=TRUE,width = 12, status = "warning", h5("No empirical PVT data loaded.")))
              ),
              conditionalPanel("output.pvt_loaded",
                fluidRow(
                  box(collapsible=TRUE,width = 12, status = "primary", solidHeader = FALSE,
                      uiOutput("pvtlap_ui_groups"),uiOutput("pvtlap_ui_subjects"),verbatimTextOutput("pvtlap_exclusions"))
                ),
                fr_fit_stat_boxes("pvtlap"),
                fr_calib_box("pvtlap"),
                fluidRow(
                  box(collapsible=TRUE,title = "Model vs Observed PVT lapses -- time series", width = 8,
                      status = "primary", solidHeader = TRUE,
                    plotlyOutput("pvtlap_plt_ts", height = "360px")),
                  box(collapsible=TRUE,title = "Predicted vs Observed", width = 4,
                      status = "warning", solidHeader = TRUE,
                    plotlyOutput("pvtlap_plt_scatter", height = "360px"))
                ),
                fluidRow(
                  box(collapsible=TRUE,title = "Residuals (Model - Observed)", width = 12,
                      status = "danger", solidHeader = TRUE,
                    plotlyOutput("pvtlap_plt_residuals", height = "200px"))
                ),
                fluidRow(
                  box(collapsible=TRUE,title = "Comparison table", width = 12,
                      status = "success", solidHeader = TRUE,
                    DTOutput("pvtlap_tbl"), br(),
                    downloadButton("pvtlap_dl", "Download CSV", class = "btn-sm btn-success"))
                )
              )
            ),

            # ---- Composite index sub-tab ----
            tabPanel("Composite index",
              conditionalPanel("!output.index_loaded",
                fluidRow(box(collapsible=TRUE,width = 12, status = "warning", h5("Load KSS (or VAS) and/or PVT observations to enable a composite index comparison.")))
              ),
              conditionalPanel("output.index_loaded",
                fluidRow(
                  box(collapsible=TRUE,width = 12, status = "primary", solidHeader = FALSE,
                      uiOutput("index_ui_groups"),uiOutput("index_ui_subjects"),verbatimTextOutput("index_exclusions"))
                ),
                fr_fit_stat_boxes("index"),
                fluidRow(
                  box(collapsible=TRUE,title = "Model vs Observed composite index -- time series", width = 8,
                      status = "primary", solidHeader = TRUE,
                    plotlyOutput("index_plt_ts", height = "360px")),
                  box(collapsible=TRUE,title = "Predicted vs Observed", width = 4,
                      status = "warning", solidHeader = TRUE,
                    plotlyOutput("index_plt_scatter", height = "360px"))
                ),
                fluidRow(
                  box(collapsible=TRUE,title = "Residuals (Model - Observed)", width = 12,
                      status = "danger", solidHeader = TRUE,
                    plotlyOutput("index_plt_residuals", height = "200px"))
                ),
                fluidRow(
                  box(collapsible=TRUE,title = "Comparison table", width = 12,
                      status = "success", solidHeader = TRUE,
                    DTOutput("index_tbl"), br(),
                    downloadButton("index_dl", "Download CSV", class = "btn-sm btn-success"))
                ),
                helpText("\"Observed\" here is not a directly measured quantity. Where a subject has both a KSS (or VAS-rescaled) and a PVT observation at the same elapsed hour, they are combined using the configured weights. Where only KSS/VAS or only PVT is loaded at all, the comparison falls back to 100% weight on that single term -- and the model series it is scored against uses that same 100% weighting, not the weights configured on the Parameters tab, so RMSE/MAE/r here may not match the mixed-weight Index shown elsewhere. Composite index itself remains author-defined and unvalidated.")
              )
            )
          )
        )
      ),  # end compare tabItem


      # ========== EXPORT TAB ==========
      tabItem(tabName = "export",
        fluidRow(
          box(collapsible=TRUE,title = "Reproducible R Script", width = 12,
              status = "primary", solidHeader = TRUE,
            downloadButton("dl_bundle", "Download full analysis (.rds)"),
            downloadButton("dl_rscript", "Download R script",
                          class = "btn-primary"),
            br(), br(),
            verbatimTextOutput("rscript_preview")
          )
        ),
        fluidRow(
          box(collapsible=TRUE,title = "Pipeline Log", width = 12,
              status = "info", solidHeader = TRUE,
            verbatimTextOutput("pipeline_log")
          )
        )
      )  # end export tabItem

    )  # end tabItems
  )  # end dashboardBody
)  # end dashboardPage



# KSS-referenced traffic-light colour for a KSS value (1 = green, 9 = red),
# with the yellow-to-orange transition centred on the 4-6 range.
fr_risk_colour <- function(kss) {
  kss <- pmin(pmax(kss, 1), 9)
  stops <- c(1, 4, 6, 9)
  cols  <- c("#2ecc71", "#f1c40f", "#e67e22", "#e74c3c")
  seg <- findInterval(kss, stops, all.inside = TRUE)
  t <- (kss - stops[seg]) / (stops[seg + 1L] - stops[seg])
  vapply(seq_along(kss), function(i)
    grDevices::colorRampPalette(cols[c(seg[i], seg[i] + 1L)])(101)[round(t[i] * 100) + 1L], character(1))
}

# Background gradient for the composite index plot. The index is on a [0,1]
# scale, so bands are coloured via its KSS-equivalent (KSS = 1 + 8*Index),
# the same transform the app already uses to map KSS onto the index scale.
fr_risk_bands <- function(y0 = 0, y1 = 1, n = 60, alpha = 0.30) {
  ys  <- seq(y0, y1, length.out = n + 1L)
  mid <- (head(ys, -1L) + tail(ys, -1L)) / 2
  cols <- fr_risk_colour(1 + 8 * mid)
  lapply(seq_len(n), function(i) list(type = "rect", layer = "below", xref = "paper",
    x0 = 0, x1 = 1, y0 = ys[i], y1 = ys[i + 1L],
    fillcolor = hex_to_rgba(cols[i], alpha), line = list(width = 0)))
}

sleep_shapes <- function(res,y0=0,y1=1) {
  s<-attr(res,"bundle")$inputs$schedule
  lapply(seq_len(nrow(s)),function(i)list(type="rect",layer="below",fillcolor="rgba(174,214,241,0.35)",
    line=list(width=0),x0=s$onset_h[i],x1=s$offset_h[i],y0=y0,y1=y1,yref="paper"))
}

server <- function(input,output,session) {
  rv<-reactiveValues(schedule_manual=fr_empty_intervals(),bundle=NULL,hash=NULL,error=NULL)
  origin<-reactive({
    if(isTRUE(input$dated)){
      req(nzchar(trimws(input$start_iso)))
      validate(need(grepl("^[0-9]{4}-[0-9]{2}-[0-9]{2}[T ][0-9]{2}:[0-9]{2}(:[0-9]{2}(\\.[0-9]+)?)?(Z|[+-][0-9]{2}:?[0-9]{2})$",trimws(input$start_iso)),
        "Study start: use ISO datetimes with an explicit UTC offset, e.g. 2026-09-16T08:00:00+02:00."))
      fr_origin("dated",start_iso=input$start_iso)
    } else {
      req(is.finite(input$t0_clock))
      validate(need(input$t0_clock>=0 && input$t0_clock<24, "Clock hour at t = 0 must be between 0 and 24."))
      fr_origin(clock=input$t0_clock)
    }
  })
  output$tbl_schedule_manual<-renderDT(datatable(rv$schedule_manual,selection="multiple",editable=TRUE,
    options=list(pageLength=15,dom="t",ordering=FALSE),rownames=FALSE,
    colnames=c("Onset (h from t=0)","Offset (h from t=0)"),class="compact stripe"))
  observeEvent(input$tbl_schedule_manual_cell_edit,{
    e<-input$tbl_schedule_manual_cell_edit
    tryCatch({d<-rv$schedule_manual;d[e$row,e$col+1]<-fr_parse_numbers(e$value,name="Sleep boundary",allow_missing=FALSE)
      rv$schedule_manual<-d},error=function(e)showNotification(conditionMessage(e),type="error"))
  })
  observeEvent(input$btn_add,{rv$schedule_manual<-rbind(rv$schedule_manual,data.frame(onset_h=input$new_onset,offset_h=input$new_offset))})
  observeEvent(input$btn_del,{sel<-input$tbl_schedule_manual_rows_selected;if(length(sel))rv$schedule_manual<-rv$schedule_manual[-sel,,drop=FALSE]})
  observeEvent(input$btn_example,{
    rv$schedule_manual<-data.frame(onset_h=15+(0:6)*24,offset_h=23+(0:6)*24)
    updateNumericInput(session,"duration",value=175);updateRadioButtons(session,"input_mode",selected="manual")
  })
  schedule_source<-reactive({
    if(input$input_mode=="manual")return(NULL)
    f<-if(input$input_mode=="file")input$schedule_file else input$actig_file
    req(!is.null(f));d<-classic_read_table(f$datapath,f$name)
    list(filename=f$name,sha256=digest::digest(file=f$datapath,algo="sha256"),raw_table=d)
  })
  active_schedule<-reactive({
    end<-input$duration
    if(input$input_mode=="manual")return(fr_intervals(rv$schedule_manual,end))
    d<-schedule_source()$raw_table
    if(input$input_mode=="file"){
      if(ncol(d)<2)stop("Schedule file needs onset and offset in its first two columns.",call.=FALSE)
      s<-data.frame(onset_h=fr_elapsed(d[[1]],input$schedule_format,origin(),name="Sleep onset"),
        offset_h=fr_elapsed(d[[2]],input$schedule_format,origin(),name="Sleep offset"))
      fr_intervals(s,end)
    } else {
      names(d)<-tolower(names(d));sc<-intersect(c("sleep_wake","state"),names(d))[1]
      if(!"datetime" %in% names(d)||is.na(sc))stop("Actigraphy needs datetime and sleep_wake (or state) columns.",call.=FALSE)
      tt<-fr_elapsed(d$datetime,"iso",origin(),name="Actigraphy time")
      fr_actigraphy(tt,d[[sc]],input$epoch_minutes/60,end)
    }
  })
  output$tbl_file_preview<-renderDT({
    x<-tryCatch(active_schedule(),error=function(e)e)
    validate(need(!inherits(x,"error"),if(inherits(x,"error"))conditionMessage(x)else""))
    datatable(x,options=list(pageLength=10,dom="tp"),rownames=FALSE,class="compact stripe")
  })
  # Read-only published values until the user deliberately enables exploration.
  parameter_box<-function(spec)renderUI({
    version<-input$model_version;req(version)
    custom<-isTRUE(input$custom_params)
    tagList(lapply(spec,function(row){
      metric<-row[1];nm<-row[2];label<-row[3];value<-mccauley_params(metric,version)[[nm]]
      desc<-unname(FR_PARAM_INFO[nm])
      note<-if(!is.na(desc))tags$div(desc,style="font-size:85%;color:#6c757d;margin:-4px 0 10px 2px;") else NULL
      if((nm%in%c("gamma","nu")&&version=="2021")||(nm=="zeta"&&version!="2024fe"))return(tagList(helpText(paste(label,": absent from this version")),note))
      if(custom){
        rng<-fr_slider_range(nm,value)
        step<-signif(max((rng[2]-rng[1])/200,1e-4),2)
        tagList(sliderInput(paste0("par_",metric,"_",nm),label,min=rng[1],max=rng[2],value=value,step=step),note)
      }
      else tagList(tags$p(tags$b(paste0(label,": ")),format(value,digits=7)),note)
    }))
  })
  output$homeostatic_controls<-parameter_box(list(c("PVT","alpha_w","PVT wake rate (1/h)"),c("PVT","alpha_s","PVT sleep rate (1/h)"),
    c("KSS","alpha_w","KSS wake rate (1/h)"),c("KSS","alpha_s","KSS sleep rate (1/h)")))
  output$allostatic_controls<-parameter_box(list(c("PVT","eta_w","Shared eta_w (1/h)"),c("PVT","beta","Shared beta"),c("PVT","zeta","Shared feedback zeta"),
    c("PVT","W_c","PVT Wc (h)"),c("KSS","W_c","KSS Wc (h)")))
  output$circadian_controls<-parameter_box(list(c("PVT","phi","Shared sine phase phi (clock h)"),c("PVT","lambda","Shared amplitude rate lambda (1/h)"),
    c("PVT","mu_s","Shared sleep offset mu_s"),c("PVT","mu_w","PVT wake offset mu_w"),c("KSS","mu_w","KSS wake offset mu_w")))
  output$inertia_controls<-parameter_box(list(c("PVT","nu","Shared inertia decay rate nu (1/h)"),c("PVT","gamma","Shared inertia scaling gamma")))
  output$kss_parameter_controls<-parameter_box(list(c("KSS","p_f","KSS floor p_f"),c("KSS","xi","KSS scaling xi")))
  output$pvt_parameter_controls<-parameter_box(list(c("PVT","p_f","PVT floor p_f"),c("PVT","xi","PVT scaling xi")))
  output$model_reference<-renderUI({
    ref<-FR_REFERENCES[[paste0("model_",input$model_version)]]
    tags$a(href=ref,target="_blank","Equations, parameters and validation in the primary paper")
  })
  output$calibration_status<-renderText(if(isTRUE(input$custom_params))"CUSTOM — UNVALIDATED parameter configuration" else "PUBLISHED DEFAULTS — calibration and validation limited to the source studies")
  overrides<-reactive({
    ans<-list(PVT=list(),KSS=list());if(!isTRUE(input$custom_params))return(ans)
    shared<-c("eta_w","beta","zeta","phi","lambda","mu_s","nu","gamma")
    specific<-c("alpha_w","alpha_s","W_c","mu_w","p_f","xi")
    for(m in names(ans))for(nm in c(shared,specific)){
      if((nm%in%c("gamma","nu")&&input$model_version=="2021")||(nm=="zeta"&&input$model_version!="2024fe"))next
      v<-input[[paste0("par_",if(nm%in%shared)"PVT"else m,"_",nm)]]
      req(!is.null(v));ans[[m]][[nm]]<-v
    }
    ans
  })
  output$weight_check<-renderText({
    tryCatch({I<-fr_index_settings(list(weights=c(input$w1,input$w2),pvt_ceiling=input$pvt_ceiling))
      sprintf("Normalized weights: KSS %.3f; PVT %.3f. Author-defined, not validated.",I$normalized_weights[1],I$normalized_weights[2])},error=function(e)conditionMessage(e))
  })
  imported<-list()
  for(m in c("KSS","VAS","PVT"))imported[[m]]<-make_wide_outcome_importer(input,output,session,NULL,rv,
    prefix=tolower(m),label=m,col_regex="",origin_fn=origin)
  raw_observations<-function(m)reactive({r<-imported[[m]]();if(is.null(r)||isTRUE(r$disabled))return(NULL)
    d<-r$observations;d[[if(m=="PVT")"PVT"else"KSS"]]<-d$value;d})
  kss_long_data<-raw_observations("KSS");vas_kss_equiv<-raw_observations("VAS");pvt_emp<-raw_observations("PVT")

  # Study duration auto-extends to cover the latest loaded measurement time +
  # 0.5h once any outcome data is loaded, so a schedule doesn't silently cut
  # off the last observation. Same auto/manual toggle pattern used elsewhere
  # (e.g. the VAS range auto-detect): stays "auto" until the user edits the
  # field themselves, so it never fights a deliberate manual value.
  duration_auto<-reactiveVal(TRUE);duration_auto_value<-reactiveVal(NULL)
  observeEvent(input$duration,{
    if(!identical(input$duration,isolate(duration_auto_value())))duration_auto(FALSE)
  },ignoreInit=TRUE)
  observe({
    times<-c()
    for(m in c("KSS","VAS","PVT")){
      r<-tryCatch(imported[[m]](),error=function(e)NULL)
      if(!is.null(r)&&!isTRUE(r$disabled)&&length(r$settings$times))times<-c(times,r$settings$times)
    }
    req(length(times)>0)
    new_dur<-max(times)+0.5
    isolate({
      if(isTRUE(duration_auto())&&!identical(input$duration,new_dur)){
        duration_auto_value(new_dur)
        updateNumericInput(session,"duration",value=new_dur)
      }
    })
  },priority=5)

  # Cross-modality agreement (GEE): unlike the CCC boxes on each comparison
  # sub-tab (which compare the MODEL to one outcome type), this compares two
  # loaded OUTCOME TYPES to each other -- e.g. does observed KSS agree with
  # observed PVT -- which only makes sense when both genuinely vary by
  # subject (see fr_inter_ccc_gee()'s header comment for why a model
  # prediction can't be used this way). Applies each metric's own current
  # group/subject picker selection (same pattern as sim_inputs_live) before
  # rescaling to a common 0-1 impairment scale, so a subject excluded on one
  # tab is naturally excluded from that tab's pairings too.
  modality_norm <- function(m) reactive({
    d <- switch(m, KSS = kss_long_data(), VAS = vas_kss_equiv(), PVT = pvt_emp())
    if (is.null(d) || !nrow(d)) return(NULL)
    key <- if (m == "PVT") "pvtlap" else tolower(m)
    groups <- sort(unique(d$group)); gv <- input[[paste0(key,"_groups_selected")]]
    sel_g <- if (length(groups) <= 1 || is.null(gv)) groups else as.character(gv)
    subjects <- sort(unique(as.character(d$subject))); sv <- input[[paste0(key,"_subjects_selected")]]
    sel_s <- if (length(subjects) <= 1 || is.null(sv)) subjects else as.character(sv)
    d <- d[d$group %in% sel_g & as.character(d$subject) %in% sel_s, , drop = FALSE]
    if (!nrow(d)) return(NULL)
    val <- switch(m,
      KSS = (d$value-1)/8,
      VAS = if (identical(input$vas_target,"01")) d$value else (d$value-1)/8,
      PVT = pmin(pmax(d$value,0)/input$pvt_ceiling,1))
    data.frame(subject = as.character(d$subject), time_h = d$time_h, value = val)
  })
  kss_norm <- modality_norm("KSS"); vas_norm <- modality_norm("VAS"); pvt_norm <- modality_norm("PVT")
  cross_modality_pairs <- list(
    list(key="kss_vas", label="KSS <-> VAS agreement", ra=kss_norm, rb=vas_norm),
    list(key="kss_pvt", label="KSS <-> PVT agreement", ra=kss_norm, rb=pvt_norm),
    list(key="vas_pvt", label="VAS <-> PVT agreement", ra=vas_norm, rb=pvt_norm))
  for (pr in cross_modality_pairs) local({
    p <- pr
    output[[paste0("vb_ccc_",p$key)]] <- renderValueBox({
      a <- p$ra(); b <- p$rb(); req(!is.null(a), !is.null(b))
      r <- fr_inter_ccc_gee(a$subject,a$time_h,a$value,b$subject,b$time_h,b$value)
      val <- r$ccc; col <- if (!is.na(val) && val < 0) "red" else "blue"
      valueBox(if (is.na(val)) "--" else sprintf("%.3f", val),
        paste0(p$label, " (", r$n_times, " occasions, ", r$n_subjects, " subjects)"),
        icon = icon("balance-scale"), color = col)
    })
  })
  output$cross_modality_section <- renderUI({
    avail <- c(KSS=!is.null(kss_long_data()), VAS=!is.null(vas_kss_equiv()), PVT=!is.null(pvt_emp()))
    if (sum(avail) < 2) return(NULL)
    present <- names(avail)[avail]
    boxes <- Filter(Negate(is.null), lapply(cross_modality_pairs, function(p) {
      pair_names <- toupper(strsplit(p$key,"_")[[1]])
      if (!all(pair_names %in% present)) return(NULL)
      column(4, valueBoxOutput(paste0("vb_ccc_",p$key), width = 12))
    }))
    fluidRow(box(collapsible=TRUE,title = "Cross-modality agreement", width = 12, status = "info", solidHeader = TRUE,
      helpText("Do the loaded outcome types agree with each other, on subjects/occasions where both were measured? This ignores the model entirely -- it's outcome vs. outcome, not outcome vs. prediction."),
      do.call(fluidRow, boxes),
      fr_stat_explainer("How is this different from the CCC boxes on the tabs below?",
        "The CCC boxes on the KSS/VAS/PVT/Index tabs compare the model's prediction to one observed outcome. This one compares two ", tags$i("observed"), " outcome types to each other -- e.g. does subjectively reported KSS sleepiness agree with objectively measured PVT lapses -- after rescaling both to the same 0-1 impairment scale (KSS: (KSS-1)/8; PVT: lapses/ceiling, clipped 0-1; VAS: per its own rescale target set above).",
        tags$br(), "Follows the GEE approach of Tsai et al. (2017, Statistics in Medicine, doi:10.1002/sim.7207) with the identity (D1) time-weighting: variance, covariance and mean-difference terms are computed within each occasion first, then summed across occasions before the ratio is taken -- so a systematic difference in typical level between occasions doesn't leak into the agreement estimate the way it would from naively pooling every (subject, occasion) pair together.",
        tags$br(), tags$b("Caveat: "), "point estimate only, no confidence interval (consistent with the other CCC boxes here). Occasions with fewer than 2 subjects providing both outcomes are silently dropped from the estimate.")
    ))
  })
  outputOptions(output, "cross_modality_section", suspendWhenHidden = FALSE)

  # Composite (Index) "observed" score. Where a subject has both a KSS (or
  # VAS-rescaled) and a PVT observation at the same elapsed hour, both terms
  # are combined with the configured weights. Where only one of the two
  # datasets is loaded at all, the composite falls back to 100% weight on
  # whichever term is available, so a comparison is still possible -- the
  # matched predicted series (Index_compare in fr_run) uses that same
  # weighting, never the globally configured one, for a fair comparison.
  # A VAS import rescaled to 0-1 (input$vas_target=="01") is not on the KSS
  # scale, so it must not be treated as a KSS-equivalent source here -- only
  # the default "kss" target (1-9) is comparable to kss_long_data().
  vas_is_kss_target<-reactive(is.null(input$vas_target)||identical(input$vas_target,"kss"))
  composite_weights<-reactive({
    has_kss<-!is.null(kss_long_data())||(vas_is_kss_target()&&!is.null(vas_kss_equiv()));has_pvt<-!is.null(pvt_emp())
    if(has_kss&&has_pvt)c(input$w1,input$w2)
    else if(has_pvt)c(0,1)
    else if(has_kss)c(1,0)
    else NULL
  })
  composite_index_settings<-reactive({
    w<-composite_weights();if(is.null(w))return(NULL)
    tryCatch(fr_index_settings(list(weights=w,pvt_ceiling=input$pvt_ceiling)),error=function(e)NULL)
  })
  composite_emp<-reactive({
    k<-kss_long_data();v<-if(vas_is_kss_target())vas_kss_equiv()else NULL;pv<-pvt_emp();kss_src<-if(!is.null(k))k else v
    I<-composite_index_settings();if(is.null(I))return(NULL)
    if(!is.null(kss_src)&&!is.null(pv)){
      m<-merge(kss_src[,c("subject","group","time_h","value")],pv[,c("subject","time_h","value")],
        by=c("subject","time_h"),suffixes=c(".kss",".pvt"))
      if(!nrow(m))return(NULL)
      val<-fr_index(m$value.pvt,m$value.kss,I)
      data.frame(subject=m$subject,group=m$group,occasion=paste0("t",m$time_h),
        time_h=m$time_h,value_raw=val,value=val)
    } else if(!is.null(pv)){
      val<-fr_index(pv$value,rep(1,nrow(pv)),I)
      data.frame(subject=pv$subject,group=pv$group,occasion=pv$occasion,time_h=pv$time_h,value_raw=val,value=val)
    } else if(!is.null(kss_src)){
      val<-fr_index(rep(0,nrow(kss_src)),kss_src$value,I)
      data.frame(subject=kss_src$subject,group=kss_src$group,occasion=kss_src$occasion,time_h=kss_src$time_h,value_raw=val,value=val)
    } else NULL
  })
  sim_inputs_live<-reactive({
    tryCatch({
      empirical<-list();all_imports<-list()
      for(m in names(imported)){
        r<-imported[[m]]();if(is.null(r))next
        all_imports[[m]]<-r
        if(isTRUE(r$disabled))next
        groups<-sort(unique(r$observations$group));key<-if(m=="PVT")"pvtlap"else tolower(m)
        gv<-input[[paste0(key,"_groups_selected")]]
        selected<-if(length(groups)<=1||is.null(gv))groups else as.character(gv)
        subjects<-sort(unique(as.character(r$observations$subject)))
        sv<-input[[paste0(key,"_subjects_selected")]]
        selected_subj<-if(length(subjects)<=1||is.null(sv))subjects else as.character(sv)
        r$selected_groups<-selected;r$selected_subjects<-selected_subj
        if(m=="VAS"&&identical(r$settings$target,"01"))r$compare_as<-"VAS01"
        empirical[[m]]<-r
      }
      comp_df<-composite_emp()
      if(!is.null(comp_df)&&nrow(comp_df)){
        groups<-sort(unique(comp_df$group))
        gv<-input$index_groups_selected
        selected<-if(length(groups)<=1||is.null(gv))groups else as.character(gv)
        subjects<-sort(unique(as.character(comp_df$subject)))
        sv<-input$index_subjects_selected
        selected_subj<-if(length(subjects)<=1||is.null(sv))subjects else as.character(sv)
        Ic<-composite_index_settings()
        all_imports$Index<-list(observations=comp_df,disabled=FALSE)
        empirical$Index<-list(observations=comp_df,selected_groups=selected,selected_subjects=selected_subj,index_settings=Ic)
      }
      p<-list(origin=origin(),end=input$duration,schedule=active_schedule(),duty=fr_empty_intervals(),
        version=input$model_version,overrides=overrides(),baseline=list(onset=input$baseline_on,offset=input$baseline_off),
        initial=list(PVT=NULL,KSS=NULL),spacing=as.numeric(input$dt_choice)/60,
        index=list(weights=c(input$w1,input$w2),pvt_ceiling=input$pvt_ceiling),
        exposure=list(metric="Index",threshold=input$index_threshold,scope="wake"),empirical=empirical,
        source=list(schedule=schedule_source(),imports=all_imports),
        assumptions=list(schedule_mode=input$input_mode,custom=isTRUE(input$custom_params),
          composite="Author-defined weights, normalization ceiling and threshold: no demonstrated validation",
          baseline="Assumed repeating pre-study sleep opportunities",display_scale=input$prediction_scale,
          vas="Endpoint rescaling only; not a validated KSS conversion",actigraphy="Experimental actual-sleep switching"))
      fr_index_settings(p$index)
      list(ok=TRUE,value=p,hash=digest::digest(p,algo="sha256"))
    },error=function(e)list(ok=FALSE,message=if(nzchar(conditionMessage(e)))conditionMessage(e)else"Complete the input controls."))
  })
  run_now<-function(x){
    rv$error<-NULL;rv$bundle<-NULL;rv$hash<-NULL
    if(!x$ok){rv$error<-x$message;return(invisible(NULL))}
    tryCatch({rv$bundle<-fr_capture_run(x$value);rv$hash<-x$hash},error=function(e){rv$error<-conditionMessage(e)})
  }
  delayed<-debounce(sim_inputs_live,500)
  observeEvent(delayed(),{if(isTRUE(input$auto_recalc))run_now(delayed())},ignoreInit=FALSE)
  observeEvent(input$btn_run,{run_now(sim_inputs_live())})
  observeEvent(input$auto_recalc,{if(isTRUE(input$auto_recalc))run_now(sim_inputs_live())},ignoreInit=TRUE)
  valid_bundle<-reactive({x<-sim_inputs_live();req(x$ok,!is.null(rv$bundle),identical(x$hash,rv$hash));rv$bundle})
  sim_results<-reactive({
    b<-valid_bundle();res<-b$results$predictions
    res$Risk<-res$Index
    if(!is.null(res$Index_compare))res$Risk_compare<-res$Index_compare
    if(identical(input$prediction_scale,"raw")){
      res$KSS<-res$KSS_raw;res$PVT<-res$PVT_raw
      if(!is.null(res$Index_compare_raw))res$Risk_compare<-res$Index_compare_raw
    }
    attr(res,"bundle")<-b;res
  })
  output$sim_done<-reactive({tryCatch({valid_bundle();TRUE},error=function(e)FALSE)})
  outputOptions(output,"sim_done",suspendWhenHidden=FALSE)
  output$sim_status_txt<-renderText({
    x<-sim_inputs_live();if(!x$ok)return(x$message)
    if(!is.null(rv$error))return(rv$error)
    if(is.null(rv$bundle)||!identical(rv$hash,x$hash))return(if(isTRUE(input$auto_recalc))"Updating…"else"Inputs changed — click Run Simulation.")
    b<-rv$bundle
    paste(if(isTRUE(input$auto_recalc))"Live"else"Manual", "—",nrow(b$results$predictions),"time points;",b$inputs$end,"hours.",
      if(isTRUE(input$custom_params))"Custom parameters: unvalidated."else"Published model defaults.",
      "Composite index: unvalidated.")
  })
  output$vb_kss<-renderValueBox({b<-valid_bundle();s<-b$results$exposure$summary
    v<-if(input$prediction_scale=="raw")s$mean_KSS_raw else s$mean_KSS_bounded
    valueBox(if(is.finite(v))round(v,2)else"--","Time-weighted mean KSS (wake)",icon=icon("moon"),color="blue")})
  output$vb_pvt<-renderValueBox({b<-valid_bundle();v<-b$results$exposure$summary$sampled_peak_PVT_raw
    if(input$prediction_scale=="bounded")v<-max(v,0)
    valueBox(if(is.finite(v))round(v,1)else"--","Sampled peak PVT lapses / 10 min (wake)",icon=icon("chart-bar"),color="blue")})
  output$vb_highrisk<-renderValueBox({b<-valid_bundle();s<-b$results$exposure$summary
    valueBox(paste(round(s$hours_above_threshold,2),"h"),paste("Wake time above illustrative index threshold",s$threshold),icon=icon("clock"),color="purple")})
  make_xaxis<-function(res)fr_xaxis_clock(res,fr_clock_step_from_input(input))
  saved_observations<-function(m){
    b<-valid_bundle();d<-b$results$comparisons[[m]]$observations;if(is.null(d))return(NULL)
    d<-d[d$exclusion=="",,drop=FALSE];d[[if(m=="PVT")"PVT"else"KSS"]]<-d$value;d
  }
  output$plt_kss_pvt<-renderPlotly({
    res<-sim_results();wake<-res;wake$KSS[!wake$wake]<-NA_real_;wake$PVT[!wake$wake]<-NA_real_
    has_pvt<-!is.null(pvt_emp())
    p<-plot_ly(wake,x=~time_h)%>%add_lines(y=~KSS,name="KSS model prediction",yaxis="y",line=list(color="#2980b9",width=2.5))
    if(has_pvt)p<-add_lines(p,y=~PVT,name="PVT model prediction",yaxis="y2",line=list(color="#c0392b",width=2,dash="dash"))
    for(m in c("KSS","VAS",if(has_pvt)"PVT")){
      key<-tolower(m)
      p<-add_outcome_overlay(p,saved_observations(m),if(m=="PVT")"PVT"else"KSS",yaxis=if(m=="PVT")"y2"else"y",
        show_indiv=isTRUE(input[[paste0(key,"_show_indiv")]]),show_mean=isTRUE(input[[paste0(key,"_show_mean")]]),
        label_prefix=if(m=="VAS")"VAS (KSS-rescaled, unvalidated)"else paste(m,"observed"),
        marker_symbol=switch(m,KSS="circle",VAS="triangle-up",PVT="square"),
        color_override=switch(m,KSS="#e74c3c",VAS="#e67e22",PVT="#8e44ad"))
    }
    lay<-list(xaxis=make_xaxis(res),yaxis=list(title="KSS",autorange=TRUE,gridcolor="#eee"),
      shapes=sleep_shapes(res),legend=list(orientation="h",y=-.30),margin=list(r=if(has_pvt)60 else 10,t=10),
      plot_bgcolor="#fefefe",paper_bgcolor="#fefefe")
    if(has_pvt)lay$yaxis2<-list(title="PVT lapses / 10 min",overlaying="y",side="right",autorange=TRUE,zeroline=FALSE)
    do.call(plotly::layout,c(list(p),lay))
  })
  output$plt_risk<-renderPlotly({
    res<-sim_results();overlay<-isTRUE(input$risk_overlay_sleep);mapped<-overlay&&identical(input$risk_overlay_scale,"risk")
    p<-plot_ly(res,x=~time_h,y=~Risk,name="Composite impairment index (unvalidated)",type="scatter",mode="lines",line=list(color="#8e44ad",width=2))
    threshold<-attr(res,"bundle")$inputs$exposure$threshold
    p<-add_lines(p,x=c(0,max(res$time_h)),y=rep(threshold,2),inherit=FALSE,name="Illustrative threshold (unvalidated)",line=list(color="#888888",dash="dot",width=1))
    if(overlay)for(m in c("KSS","VAS")){
      d<-saved_observations(m);if(is.null(d))next
      d$plotval<-if(mapped)(d$KSS-1)/8 else d$KSS;key<-tolower(m)
      p<-add_outcome_overlay(p,d,"plotval",yaxis=if(mapped)"y"else"y2",show_indiv=isTRUE(input[[paste0(key,"_show_indiv")]]),
        show_mean=isTRUE(input[[paste0(key,"_show_mean")]]),label_prefix=if(m=="VAS")"VAS (KSS-rescaled, unvalidated)"else"KSS observed",
        marker_symbol=if(m=="VAS")"triangle-up"else"circle",color_override=if(m=="VAS")"#e67e22"else"#e74c3c")
    }
    if(isTRUE(input$risk_overlay_pvt)){
      d<-saved_observations("PVT")
      pvt_ceiling<-attr(res,"bundle")$inputs$index$pvt_ceiling
      if(!is.null(d)&&!is.null(pvt_ceiling)){
        d$plotval<-pmin(pmax(d$PVT,0)/pvt_ceiling,1)
        p<-add_outcome_overlay(p,d,"plotval",yaxis="y",show_indiv=isTRUE(input$pvt_show_indiv),
          show_mean=isTRUE(input$pvt_show_mean),label_prefix="PVT observed (mapped, unvalidated)",marker_symbol="square",
          color_override="#8e44ad")
      }
    }
    lay<-list(xaxis=make_xaxis(res),yaxis=list(title="Index [0–1], not probability",range=c(0,1.02)),shapes=c(fr_risk_bands(),sleep_shapes(res)),
      legend=list(orientation="h",y=-.35),margin=list(t=10,r=if(overlay&&!mapped)60 else 10),plot_bgcolor="#fefefe",paper_bgcolor="#fefefe")
    if(overlay&&!mapped)lay$yaxis2<-list(title="Observed sleepiness",overlaying="y",side="right",autorange=TRUE)
    do.call(plotly::layout,c(list(p),lay))
  })
  output$plt_states<-renderPlotly({
    res<-sim_results();req(input$state_metric,input$state_variable)
    col<-if(input$state_variable=="c")"c"else paste(input$state_variable,input$state_metric,sep="_")
    plot_ly(res,x=~time_h,y=res[[col]],type="scatter",mode="lines",name=col,line=list(color="#2980b9",width=2))%>%
      layout(xaxis=make_xaxis(res),yaxis=list(title=col,autorange=TRUE),shapes=sleep_shapes(res),margin=list(t=10),plot_bgcolor="#fefefe",paper_bgcolor="#fefefe")
  })
  output$tbl_preds<-renderDT({res<-sim_results();d<-res[res$wake&abs(res$time_h-round(res$time_h))<1e-9,c("time_h","clock_h","KSS","PVT","Risk")]
    names(d)<-c("Time (h)","Clock (h)","KSS","PVT lapses","Index (unvalidated)")
    datatable(round(d,3),options=list(pageLength=20,scrollX=TRUE,dom="tp"),rownames=FALSE,class="compact stripe")})
  output$dl_csv<-downloadHandler(filename=function()paste0("farima_predictions_",Sys.Date(),".csv"),content=function(file){write.csv(valid_bundle()$results$predictions,file,row.names=FALSE)})
  make_comparison_outputs(input,output,session,"kss","KSS",kss_long_data,sim_results,"KSS","KSS",c(1,9))
  make_comparison_outputs(input,output,session,"vas","VAS / other-scale overlay (unvalidated)",vas_kss_equiv,sim_results,"KSS","KSS",c(1,9))
  make_comparison_outputs(input,output,session,"pvtlap","PVT lapses",pvt_emp,sim_results,"PVT","PVT",c(0,30))
  # VAS is scored against the KSS prediction columns (fr_compare has no "VAS"
  # case and falls through to its KSS default), so its calibration search
  # also targets the KSS arm's parameters, not a separate VAS set.
  make_calibration_suggester(input,output,session,"kss","KSS","KSS",FR_CALIB_PARAM_CHOICES,sim_inputs_live)
  make_calibration_suggester(input,output,session,"vas","KSS","VAS",FR_CALIB_PARAM_CHOICES,sim_inputs_live)
  make_calibration_suggester(input,output,session,"pvtlap","PVT","PVT",FR_CALIB_PARAM_CHOICES,sim_inputs_live)
  output$index_loaded<-reactive({d<-composite_emp();!is.null(d)&&nrow(d)>0})
  outputOptions(output,"index_loaded",suspendWhenHidden=FALSE)
  make_comparison_outputs(input,output,session,"index","Composite index (unvalidated)",composite_emp,sim_results,"Risk_compare","Risk_compare",c(0,1))
  output$dl_rscript<-downloadHandler(filename=function()paste0("farima_script_",Sys.Date(),".R"),content=function(file)writeLines(fr_export_script(valid_bundle()),file))
  output$dl_bundle<-downloadHandler(filename=function()paste0("farima_analysis_",Sys.Date(),".rds"),content=function(file)saveRDS(valid_bundle(),file))
  output$rscript_preview<-renderText({b<-valid_bundle();paste(head(strsplit(fr_export_script(b),"\n")[[1]],45),collapse="\n")})
  output$pipeline_log<-renderText({b<-valid_bundle();paste("Completed:",b$metadata$created_utc,
    "\nModel:",b$inputs$version,"\nCode SHA-256:",b$metadata$code_sha256,
    "\nSleep opportunities:",nrow(b$inputs$schedule),"\nDuration:",b$inputs$end,"hours",
    "\nParameter configuration:",if(b$inputs$assumptions$custom)"CUSTOM — unvalidated"else"Published defaults",
    "\nComposite index: author-defined, not validated.",
    "\nExport contains complete functions, parameters, source tables, mapped observations, comparison groups and reproduction checks.")})
}

app<-shinyApp(ui=ui,server=server)
if(isTRUE(getOption("farimato.launch",TRUE)))app
