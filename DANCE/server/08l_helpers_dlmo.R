# ==============================================================================
# server/08l_helpers_dlmo.R -- dim-light melatonin onset (DLMO)
#
# Pure functions, no Shiny: the module (server/76_dlmo.R), the exported script
# and the tests all call these. Nothing here knows about any particular data
# set: a profile is a vector of sampling times (hours) and concentrations.
#
# METHODS AND THEIR SOURCES
#   fixed threshold    the interpolated time at which the concentration first
#                      rises above a set threshold and stays above it (Lewy &
#                      Sack 1989; Benloucif et al. 2008, J Biol Rhythms 23:3).
#                      Linear interpolation between the two samples that
#                      bracket the crossing. In-house: no package does this.
#   variable threshold the same crossing rule, the threshold being the mean of
#                      the first n baseline samples + k SD (Voultsios, Kennaway
#                      & Dawson 1997, J Biol Rhythms 12:457; n = 3, k = 2).
#   hockey-stick       Danilenko et al. 2014, Chronobiol Int 31:349, through
#                      the dlmoR package (Thalji & Spitschan 2026, J Biol
#                      Rhythms 41:301). NOT reimplemented here: when dlmoR is
#                      not installed the method is unavailable, and says so.
# CALCULATIONS_LOG.md section 23 states every rule and why.
# ==============================================================================

DANCE_DLMO_METHOD_LABEL <- c(fixed = "Fixed threshold", variable = "Variable threshold",
                             hockey = "Hockey-stick (dlmoR)")

# Statuses. Only "determined" enters the comparisons; every other one is
# counted and reported, never given a time it does not have.
DANCE_DLMO_STATUS_LABEL <- c(
  determined      = "determined",
  before_window   = "above threshold from the first sample (onset before the window)",
  no_onset        = "no sustained rise above threshold (onset after the window, or none)",
  unconfirmed     = "crossing too close to the last sample to confirm",
  gap             = "crossing between samples further apart than the interpolation limit",
  no_baseline     = "too few baseline samples for the variable threshold",
  few_samples     = "fewer than two observed samples",
  unavailable     = "dlmoR is not installed",
  failed          = "the hockey-stick fit failed")

# ---- PROFILES ------------------------------------------------------------------
# Each row of Y is one recording; its columns are sampled at `times` (hours,
# increasing). A row whose scheduled times contain a gap longer than gap_h is
# several profiles (two nights in one row, for instance), numbered in order. A
# night column, when given, labels the rows instead (one row per night).
# Returns one list per profile: id, night, group, row, t, y (missing kept as NA).
dance_dlmo_profiles <- function(Y, times, ids, night = NULL, group = NULL, gap_h = 6) {
  Y <- as.matrix(Y); times <- as.numeric(times)
  stopifnot(ncol(Y) == length(times))
  o <- order(times); times <- times[o]; Y <- Y[, o, drop = FALSE]
  seg <- cumsum(c(1L, diff(times) > gap_h))
  nseg <- max(seg)
  out <- list()
  for (i in seq_len(nrow(Y))) for (s in seq_len(nseg)) {
    k <- seg == s
    nl <- c(if (!is.null(night)) as.character(night[i]), if (nseg > 1L) sprintf("Night %d", s))
    out[[length(out) + 1L]] <- list(
      id = as.character(ids[i]), row = i, segment = s,
      night = if (length(nl)) paste(nl, collapse = " / ") else "1",
      group = if (!is.null(group)) as.character(group[i]) else NA_character_,
      t = times[k], y = as.numeric(Y[i, k]))
  }
  out
}

# A common evening axis: each profile is moved by whole days so its first
# sample falls in [anchor, anchor + 24). Times within a profile keep their
# spacing, so a second night recorded 24 h after the first is read on the same
# clock, and 01:00 is 25, after 23:00, rather than 1, before it.
dance_dlmo_day_shift <- function(t_first, anchor) 24 * floor((t_first - anchor) / 24)

# ---- VALUES BELOW THE LIMIT OF QUANTIFICATION ------------------------------------
# When an LLOQ is given, an observed value below it is replaced by LLOQ/2, the
# LLOQ, or 0 (rule "half", "lloq", "zero"). Without an LLOQ, values are used as
# entered.
dance_dlmo_lloq <- function(y, lloq = NA_real_, rule = c("half", "lloq", "zero")) {
  rule <- match.arg(rule)
  if (!is.finite(lloq)) return(list(y = y, n_below = 0L))
  b <- is.finite(y) & y < lloq
  y[b] <- switch(rule, half = lloq / 2, lloq = lloq, zero = 0)
  list(y = y, n_below = sum(b))
}

# ---- THE THRESHOLD CROSSING ------------------------------------------------------
# The first sample above the threshold (strictly) whose `confirm` following
# samples are above it too, preceded by a sample at or below it. A rise that
# falls back sooner is skipped (counted in `spikes`). Onset is interpolated
# linearly between the two observed samples that bracket the crossing;
# missing samples are never filled in.
dance_dlmo_cross <- function(t, y, threshold, confirm = 1L, max_gap = Inf) {
  res <- function(status, time = NA_real_, i0 = NA_integer_, i1 = NA_integer_, spikes = 0L)
    list(time = time, status = status, threshold = threshold, t_before = if (is.na(i0)) NA_real_ else t[i0],
         t_after = if (is.na(i1)) NA_real_ else t[i1], spikes = spikes)
  ok <- is.finite(t) & is.finite(y)
  t <- t[ok]; y <- y[ok]
  n <- length(y)
  if (n < 2L || !is.finite(threshold)) return(res("few_samples"))
  above <- y > threshold
  spikes <- 0L
  for (i in seq_len(n)) {
    if (!above[i] || (i > 1L && above[i - 1L])) next
    later <- if (confirm > 0L && i < n) (i + 1L):min(n, i + confirm) else integer(0)
    if (length(later) && !all(above[later])) { spikes <- spikes + 1L; next }
    if (i == 1L) return(res("before_window", i1 = 1L, spikes = spikes))
    tc <- t[i - 1L] + (threshold - y[i - 1L]) * (t[i] - t[i - 1L]) / (y[i] - y[i - 1L])
    st <- if (t[i] - t[i - 1L] > max_gap) "gap"
          else if (length(later) < confirm) "unconfirmed" else "determined"
    return(res(st, tc, i - 1L, i, spikes))
  }
  res("no_onset", spikes = spikes)
}

# The variable threshold: mean + k SD of the first n_base OBSERVED samples.
dance_dlmo_variable_threshold <- function(y, n_base = 3L, k = 2) {
  yb <- y[is.finite(y)]
  if (length(yb) < max(2L, n_base)) return(NA_real_)
  yb <- yb[seq_len(n_base)]
  mean(yb) + k * stats::sd(yb)
}

# ---- HOCKEY-STICK, THROUGH dlmoR ------------------------------------------------
# dlmoR works on date-times; the hours are placed on an arbitrary reference
# date and read back against it. A result that comes back as a time of day only
# is put on the day that lies within the profile's sampled window.
dance_dlmo_as_hours <- function(x, origin, t_range) {
  if (is.null(x) || !length(x)) return(NA_real_)
  x <- x[1]
  h <- if (inherits(x, "POSIXt")) as.numeric(difftime(x, origin, units = "hours"))
       else if (inherits(x, "difftime")) as.numeric(x, units = "hours")
       else if (is.character(x)) {
         p <- suppressWarnings(as.numeric(strsplit(trimws(x), ":")[[1]]))
         if (length(p) < 2L || any(!is.finite(p))) NA_real_ else sum(p * c(1, 1 / 60, 1 / 3600)[seq_along(p)])
       } else suppressWarnings(as.numeric(x))
  if (!is.finite(h)) return(NA_real_)
  if (h >= t_range[1] - 12 && h <= t_range[2] + 12) return(h)
  # a time of day: the day that puts it nearest the window
  cand <- (h %% 24) + 24 * (floor(t_range[1] / 24) + (-1:2))
  dist <- pmax(t_range[1] - cand, cand - t_range[2], 0)
  cand[which.min(dist)]
}

dance_dlmo_hockey <- function(t, y, threshold, interval_limit = 2, fine = TRUE) {
  out <- list(time = NA_real_, status = "unavailable", threshold = threshold, note = "")
  if (!requireNamespace("dlmoR", quietly = TRUE)) return(out)
  ok <- is.finite(t) & is.finite(y)
  if (sum(ok) < 3L) { out$status <- "few_samples"; return(out) }
  origin <- as.POSIXct("2000-01-01 00:00:00", tz = "UTC")
  df <- data.frame(datetime = origin + t[ok] * 3600, melatonin = y[ok])
  msgs <- character(0)
  r <- tryCatch(withCallingHandlers(
    dlmoR::calculate_dlmo(data = df, threshold = threshold, interval_limit = interval_limit,
                          fine_flag = fine),
    message = function(m) { msgs <<- c(msgs, trimws(conditionMessage(m))); invokeRestart("muffleMessage") },
    warning = function(w) { msgs <<- c(msgs, trimws(conditionMessage(w))); invokeRestart("muffleWarning") }),
    error = function(e) e)
  if (inherits(r, "error")) {
    out$status <- "failed"; out$note <- conditionMessage(r); return(out)
  }
  tm <- if (fine) r$dlmo$fine$time %||% r$dlmo$coarse$time else r$dlmo$coarse$time
  tm <- tm %||% r$dlmo_time
  h <- dance_dlmo_as_hours(tm, origin, range(t[ok]))
  out$time <- h
  out$status <- if (is.finite(h)) "determined" else "failed"
  out$note <- paste(unique(msgs), collapse = "; ")
  out$fit <- r
  out
}

# ---- ONE PROFILE, ONE METHOD -------------------------------------------------------
dance_dlmo_estimate <- function(p, method, opt) {
  lq <- dance_dlmo_lloq(p$y, opt$lloq %||% NA_real_, opt$lloq_rule %||% "half")
  y <- lq$y
  thr <- switch(method,
    fixed = opt$threshold,
    variable = dance_dlmo_variable_threshold(y, opt$n_base %||% 3L, opt$k_sd %||% 2),
    hockey = opt$threshold)
  r <- if (method == "variable" && !is.finite(thr)) list(time = NA_real_, status = "no_baseline", threshold = NA_real_)
       else if (method == "hockey") dance_dlmo_hockey(p$t, y, thr, opt$interval_limit %||% 2, isTRUE(opt$fine %||% TRUE))
       else dance_dlmo_cross(p$t, y, thr, opt$confirm %||% 1L, opt$max_gap %||% Inf)
  data.frame(id = p$id, night = p$night, group = p$group, method = method,
             threshold = r$threshold %||% thr, dlmo_h = r$time, status = r$status,
             n_samples = sum(is.finite(p$t) & is.finite(p$y)), n_below_lloq = lq$n_below,
             spikes_skipped = r$spikes %||% NA_integer_,
             t_before = r$t_before %||% NA_real_, t_after = r$t_after %||% NA_real_,
             note = r$note %||% "", stringsAsFactors = FALSE)
}

# ---- HABITUAL BEDTIME AND THE PHASE ANGLE -------------------------------------------
# Bedtime as a time of day in hours [0, 24): "23:30", "23.5", 23.5, an Excel day
# fraction (0.979), or a date-time.
dance_dlmo_clock_hours <- function(x) {
  if (inherits(x, "POSIXt")) {
    lt <- as.POSIXlt(x)
    return(lt$hour + lt$min / 60 + lt$sec / 3600)
  }
  if (inherits(x, "difftime")) return(as.numeric(x, units = "hours") %% 24)
  x <- as.character(x)
  vapply(x, function(s) {
    s <- trimws(s)
    if (is.na(s) || !nzchar(s)) return(NA_real_)
    if (grepl(":", s, fixed = TRUE)) {
      p <- suppressWarnings(as.numeric(strsplit(s, ":", fixed = TRUE)[[1]]))
      if (length(p) < 2L || any(!is.finite(p))) return(NA_real_)
      h <- sum(p * c(1, 1 / 60, 1 / 3600)[seq_along(p)])
    } else {
      h <- suppressWarnings(as.numeric(s))
      if (!is.finite(h)) return(NA_real_)
      if (h > 0 && h < 1 && !grepl("^[0-9]+$", s)) h <- h * 24   # a day fraction
    }
    if (h < 0 || h > 24) NA_real_ else h %% 24
  }, numeric(1), USE.NAMES = FALSE)
}

# Phase angle of entrainment = bedtime - DLMO, in hours: positive when the
# onset precedes bedtime. The bedtime is the one within 12 h of the onset.
dance_dlmo_phase_angle <- function(dlmo_h, bed_clock) {
  bed <- bed_clock + 24 * round((dlmo_h - bed_clock) / 24)
  bed - dlmo_h
}

# ---- GROUPS AND NIGHTS ------------------------------------------------------------
# Profiles without a determined onset are excluded and counted per cell. With a
# within-participant factor (nights) the comparison is paired: a paired t test
# for two levels and no group factor, otherwise a linear mixed model with a
# random intercept per participant (lmerTest, Satterthwaite df). Between groups
# only: Welch's t test (two groups) or Welch's one-way ANOVA.
dance_dlmo_describe <- function(d, value, by) {
  key <- if (length(by)) interaction(d[by], drop = TRUE, sep = " / ") else factor(rep("All", nrow(d)))
  do.call(rbind, lapply(levels(key), function(k) {
    s <- d[key == k, , drop = FALSE]
    v <- s[[value]][s$status == "determined" & is.finite(s[[value]])]
    data.frame(cell = k, n_profiles = nrow(s), n_determined = length(v),
               n_before = sum(s$status == "before_window"), n_after = sum(s$status == "no_onset"),
               n_other = nrow(s) - length(v) - sum(s$status %in% c("before_window", "no_onset")),
               mean = if (length(v)) mean(v) else NA_real_,
               sd = if (length(v) > 1L) stats::sd(v) else NA_real_,
               median = if (length(v)) stats::median(v) else NA_real_,
               min = if (length(v)) min(v) else NA_real_, max = if (length(v)) max(v) else NA_real_,
               stringsAsFactors = FALSE)
  }))
}

dance_dlmo_compare <- function(d, value = "dlmo_h", between = FALSE, within = FALSE, conf = 0.95) {
  d <- d[d$status == "determined" & is.finite(d[[value]]), , drop = FALSE]
  d$v <- d[[value]]
  none <- function(msg) list(ok = FALSE, message = msg, tests = NULL)
  if (!between && !within) return(none("No factor to compare."))
  if (nrow(d) < 3L) return(none("Fewer than three determined onsets."))
  d$group <- factor(d$group); d$night <- factor(d$night)
  if (between && nlevels(d$group) < 2L) return(none("Only one group has determined onsets."))
  if (within && nlevels(d$night) < 2L) return(none("Only one night has determined onsets."))
  if (!within && any(duplicated(d$id)))
    return(none("Participants contribute several profiles; compare the nights as a within-participant factor."))
  row <- function(effect, test, est = NA, lo = NA, hi = NA, stat = NA, df1 = NA, df2 = NA, p = NA, n = NA)
    data.frame(effect = effect, test = test, estimate = est, lo = lo, hi = hi, statistic = stat,
               df1 = df1, df2 = df2, p = p, n = n, stringsAsFactors = FALSE)
  if (within && !between && nlevels(d$night) == 2L) {
    w <- stats::reshape(d[, c("id", "night", "v")], idvar = "id", timevar = "night", direction = "wide")
    lv <- levels(d$night); a <- w[[paste0("v.", lv[2])]]; b <- w[[paste0("v.", lv[1])]]
    cc <- is.finite(a) & is.finite(b)
    if (sum(cc) < 2L) return(none("Fewer than two participants have an onset on both nights."))
    tt <- stats::t.test(a[cc], b[cc], paired = TRUE, conf.level = conf)
    return(list(ok = TRUE, method = "paired t test", n_pairs = sum(cc), n_unpaired = length(unique(d$id)) - sum(cc),
                tests = row(sprintf("%s - %s", lv[2], lv[1]), "paired t", unname(tt$estimate),
                            tt$conf.int[1], tt$conf.int[2], unname(tt$statistic), unname(tt$parameter),
                            NA, tt$p.value, sum(cc))))
  }
  if (!within) {
    if (nlevels(d$group) == 2L) {
      lv <- levels(d$group)
      tt <- stats::t.test(d$v[d$group == lv[2]], d$v[d$group == lv[1]], conf.level = conf)
      return(list(ok = TRUE, method = "Welch t test",
                  tests = row(sprintf("%s - %s", lv[2], lv[1]), "Welch t", unname(diff(rev(tt$estimate))),
                              tt$conf.int[1], tt$conf.int[2], unname(tt$statistic), unname(tt$parameter),
                              NA, tt$p.value, nrow(d))))
    }
    ow <- stats::oneway.test(v ~ group, data = d)
    return(list(ok = TRUE, method = "Welch one-way ANOVA",
                tests = row("group", "Welch F", stat = unname(ow$statistic), df1 = unname(ow$parameter[1]),
                            df2 = unname(ow$parameter[2]), p = ow$p.value, n = nrow(d))))
  }
  if (!requireNamespace("lmerTest", quietly = TRUE)) return(none("The mixed model needs the lmerTest package."))
  fml <- if (between) v ~ group * night + (1 | id) else v ~ night + (1 | id)
  m <- tryCatch(suppressMessages(lmerTest::lmer(fml, data = d)), error = function(e) e)
  if (inherits(m, "error")) return(none(paste("The mixed model failed:", conditionMessage(m))))
  a <- as.data.frame(stats::anova(m, ddf = "Satterthwaite"))
  tests <- do.call(rbind, lapply(seq_len(nrow(a)), function(i)
    row(rownames(a)[i], "F (Satterthwaite)", stat = a[i, "F value"], df1 = a[i, "NumDF"],
        df2 = a[i, "DenDF"], p = a[i, "Pr(>F)"], n = nrow(d))))
  # the night difference within each group (or overall), from the model
  em <- NULL
  if (nlevels(d$night) == 2L) {
    lv <- levels(d$night)
    gl <- if (between) levels(d$group) else NA
    em <- do.call(rbind, lapply(gl, function(g) {
      nd <- data.frame(night = factor(lv, levels = lv), group = if (between) factor(g, levels = levels(d$group)) else NA)
      X <- stats::model.matrix(stats::delete.response(stats::terms(lme4::nobars(fml))), nd)
      L <- X[2, names(lme4::fixef(m))] - X[1, names(lme4::fixef(m))]
      ct <- lmerTest::contest(m, matrix(L, nrow = 1), joint = FALSE, confint = TRUE, level = conf)
      row(sprintf("%s - %s%s", lv[2], lv[1], if (between) paste0(" in ", g) else ""), "contrast (Satterthwaite)",
          ct$Estimate, ct$lower, ct$upper, ct[["t value"]], ct$df, NA, ct[["Pr(>|t|)"]],
          sum(if (between) d$group == g else TRUE))
    }))
  }
  list(ok = TRUE, method = "linear mixed model, random intercept per participant",
       tests = rbind(tests, em), model = m)
}

# Two methods on the same profiles: the mean and SD of their difference.
dance_dlmo_agreement <- function(tab) {
  ms <- unique(tab$method)
  if (length(ms) < 2L) return(NULL)
  key <- paste(tab$id, tab$night, sep = "\r")
  prs <- utils::combn(ms, 2, simplify = FALSE)
  do.call(rbind, lapply(prs, function(pr) {
    a <- tab[tab$method == pr[2], ]; b <- tab[tab$method == pr[1], ]
    ka <- paste(a$id, a$night, sep = "\r"); kb <- paste(b$id, b$night, sep = "\r")
    k <- intersect(ka[a$status == "determined"], kb[b$status == "determined"])
    dd <- a$dlmo_h[match(k, ka)] - b$dlmo_h[match(k, kb)]
    data.frame(comparison = sprintf("%s - %s", DANCE_DLMO_METHOD_LABEL[pr[2]], DANCE_DLMO_METHOD_LABEL[pr[1]]),
               n_both = length(dd), mean_diff_h = if (length(dd)) mean(dd) else NA_real_,
               sd_diff_h = if (length(dd) > 1L) stats::sd(dd) else NA_real_,
               max_abs_diff_h = if (length(dd)) max(abs(dd)) else NA_real_, stringsAsFactors = FALSE)
  }))
}

# ---- THE WHOLE RUN ----------------------------------------------------------------
# profiles: from dance_dlmo_profiles(); methods: any of fixed, variable, hockey;
# opt: threshold, confirm, max_gap, n_base, k_sd, lloq, lloq_rule,
# interval_limit, fine, anchor; bedtime: NULL, or a time of day per row of the
# original data (NA where unknown).
dance_dlmo_run <- function(profiles, methods, opt, bedtime = NULL) {
  anchor <- opt$anchor %||% floor(min(vapply(profiles, function(p) min(p$t) %% 24, numeric(1))))
  tab <- do.call(rbind, lapply(profiles, function(p) {
    sh <- dance_dlmo_day_shift(min(p$t), anchor)
    p$t <- p$t - sh
    r <- do.call(rbind, lapply(methods, function(m) dance_dlmo_estimate(p, m, opt)))
    r$bedtime_h <- if (!is.null(bedtime)) bedtime[p$row] else NA_real_
    r$day_shift <- sh
    r
  }))
  tab$dlmo_clock <- dance_clock_label(tab$dlmo_h, 24, show_day = FALSE)
  tab$phase_angle_h <- if (!is.null(bedtime)) dance_dlmo_phase_angle(tab$dlmo_h, tab$bedtime_h) else NA_real_
  rownames(tab) <- NULL
  attr(tab, "anchor") <- anchor
  tab
}

# ---- REPORTING ---------------------------------------------------------------------
# r: the run as the module stores it (tab, methods, opt, has_night, has_group,
# has_bed, time_src, units, gap_h, bed_src, ...). Every rule the run applied,
# in words, for the summary, the report and the script.
dance_dlmo_settings_text <- function(r) {
  o <- r$opt; u <- r$units %||% "pg/mL"
  c(sprintf("Profiles: %d, from %d recording(s); a gap of more than %s h between scheduled samples starts a new profile%s.",
            length(r$profiles), r$n_rows, format(r$gap_h), if (r$has_night_var) sprintf("; nights named by '%s'", r$night_var) else ""),
    sprintf("Sampling times: %s; read on one evening clock from %s (each profile moved by whole days).",
            if (identical(r$time_src, "_manual_")) "entered by hand" else "clock times parsed from the column names",
            dance_clock_label(r$anchor, 24, show_day = FALSE)),
    if ("fixed" %in% r$methods) sprintf("Fixed threshold: %s %s.", format(o$threshold), u),
    if ("variable" %in% r$methods) sprintf("Variable threshold: mean + %s SD of the first %d observed samples of each profile.",
                                            format(o$k_sd), as.integer(o$n_base)),
    if (any(c("fixed", "variable") %in% r$methods))
      sprintf("Crossing: the first sample above the threshold with %d following sample(s) also above it; onset interpolated linearly between the two observed samples bracketing the crossing, if they are at most %s h apart. Missing samples are not filled in.",
              as.integer(o$confirm), format(o$max_gap)),
    if ("hockey" %in% r$methods) sprintf("Hockey-stick: dlmoR %s (Thalji & Spitschan, 2026), threshold %s %s, repeated crossings within %s h treated as one rise, fine grid search.",
                                          r$dlmoR_version %||% "(not installed)", format(o$threshold), u, format(o$interval_limit)),
    if (is.finite(o$lloq %||% NA_real_)) sprintf("Values below the LLOQ (%s %s) set to %s.", format(o$lloq), u,
                                                 switch(o$lloq_rule, half = "LLOQ/2", lloq = "the LLOQ", zero = "0"))
    else "No limit of quantification given: values used as entered.",
    if (r$has_bed) sprintf("Habitual bedtime: %s. Phase angle = bedtime - DLMO (positive: onset before bedtime), using the bedtime within 12 h of the onset.",
                           if (identical(r$bed_src, "value")) r$bed_value else sprintf("column '%s'", r$bed_var))
    else "Habitual bedtime not given: no phase angle.",
    "Profiles without a determined onset are excluded from the comparisons and counted per cell.")
}

# Descriptives and tests for one method's onsets (and phase angles).
dance_dlmo_comparisons <- function(r, method, conf = 0.95) {
  d <- r$tab[r$tab$method == method, , drop = FALSE]
  by <- c(if (r$has_group) "group", if (r$has_night) "night")
  out <- list(method = method,
              desc = dance_dlmo_describe(d, "dlmo_h", by),
              test = dance_dlmo_compare(d, "dlmo_h", r$has_group, r$has_night, conf))
  if (r$has_bed) {
    out$pa_desc <- dance_dlmo_describe(d, "phase_angle_h", by)
    out$pa_test <- dance_dlmo_compare(d, "phase_angle_h", r$has_group, r$has_night, conf)
  }
  out
}

# A difference in hours as "+1.25 h (75 min)".
dance_dlmo_fmt_diff <- function(h, digits = 2)
  ifelse(is.finite(h), sprintf("%+.*f h (%+.0f min)", digits, h, 60 * h), "NA")

dance_dlmo_test_text <- function(tr) {
  p <- if (is.finite(tr$p)) (if (tr$p < .001) "p < .001" else sprintf("p = %s", sub("^0", "", sprintf("%.3f", tr$p)))) else ""
  st <- if (tr$test %in% c("paired t", "Welch t", "contrast (Satterthwaite)"))
          sprintf("t(%.1f) = %.2f", tr$df1, tr$statistic)
        else sprintf("F(%s, %.1f) = %.2f", format(round(tr$df1, 1)), tr$df2, tr$statistic)
  paste(c(st, p), collapse = ", ")
}

# Markdown for the APA report: Methods and Results paragraphs.
dance_dlmo_apa <- function(r) {
  L <- character(0); add <- function(...) L <<- c(L, paste0(...))
  ml <- DANCE_DLMO_METHOD_LABEL[r$methods]
  cite <- c(fixed = "the fixed-threshold method (Lewy & Sack, 1989; Benloucif et al., 2008)",
            variable = "the variable-threshold method (Voultsios et al., 1997)",
            hockey = "the hockey-stick method (Danilenko et al., 2014) as implemented in dlmoR (Thalji & Spitschan, 2026)")
  add(sprintf("DLMO was estimated in %d profile(s) by %s. ", length(r$profiles), paste(cite[r$methods], collapse = "; ")),
      paste(dance_dlmo_settings_text(r), collapse = " "))
  add("")
  for (m in r$methods) {
    cp <- dance_dlmo_comparisons(r, m)
    ds <- cp$desc
    cells <- vapply(seq_len(nrow(ds)), function(i) sprintf("%s: %s (SD %s h), %d of %d profiles%s",
      ds$cell[i], if (is.finite(ds$mean[i])) dance_clock_label(ds$mean[i], 24, show_day = FALSE) else "NA",
      if (is.finite(ds$sd[i])) sprintf("%.2f", ds$sd[i]) else "NA", ds$n_determined[i], ds$n_profiles[i],
      if (ds$n_before[i] + ds$n_after[i] + ds$n_other[i] > 0)
        sprintf("; not determined: %d before the window, %d no onset, %d other", ds$n_before[i], ds$n_after[i], ds$n_other[i])
      else ""), "")
    add(sprintf("*%s.* Mean DLMO %s.", ml[[m]], paste(cells, collapse = "; ")))
    if (isTRUE(cp$test$ok)) {
      tt <- cp$test$tests
      add(sprintf("Comparison (%s%s): %s.", cp$test$method,
                  if (!is.null(cp$test$n_pairs)) sprintf(", %d participants with an onset on both nights, %d without", cp$test$n_pairs, cp$test$n_unpaired) else "",
                  paste(vapply(seq_len(nrow(tt)), function(i) sprintf("%s%s, %s", tt$effect[i],
                    if (is.finite(tt$estimate[i])) sprintf(" %s [%s, %s]", dance_dlmo_fmt_diff(tt$estimate[i]),
                                                            sprintf("%.2f", tt$lo[i]), sprintf("%.2f", tt$hi[i])) else "",
                    dance_dlmo_test_text(tt[i, ])), ""), collapse = "; ")))
    } else if (r$has_group || r$has_night) add(sprintf("No comparison: %s", cp$test$message))
    if (r$has_bed) {
      pd <- cp$pa_desc
      add(sprintf("Phase angle (bedtime - DLMO): %s.", paste(vapply(seq_len(nrow(pd)), function(i)
        sprintf("%s: M = %s h (SD %s), n = %d", pd$cell[i], sprintf("%.2f", pd$mean[i]),
                if (is.finite(pd$sd[i])) sprintf("%.2f", pd$sd[i]) else "NA", pd$n_determined[i]), ""), collapse = "; ")))
      if (isTRUE(cp$pa_test$ok)) {
        tt <- cp$pa_test$tests
        add(sprintf("Phase angle comparison (%s): %s.", cp$pa_test$method,
                    paste(vapply(seq_len(nrow(tt)), function(i) sprintf("%s, %s", tt$effect[i], dance_dlmo_test_text(tt[i, ])), ""),
                          collapse = "; ")))
      }
    }
    add("")
  }
  ag <- dance_dlmo_agreement(r$tab)
  if (!is.null(ag)) for (i in seq_len(nrow(ag)))
    add(sprintf("Agreement, %s: mean difference %s (SD %.2f h) over %d profiles with both onsets.",
                ag$comparison[i], dance_dlmo_fmt_diff(ag$mean_diff_h[i]), ag$sd_diff_h[i], ag$n_both[i]))
  L
}
