# ==============================================================================
# server/06_helpers_mixed.R — mixed (between x within) designs
# ==============================================================================
# A design with one BETWEEN-subject factor and one WITHIN-subject factor cannot
# be analysed by either of the app's existing kernels. perform_functional_anova()
# is one-way between (whole curves permuted across groups); perform_rm_fanova()
# is one-factor repeated measures (condition labels permuted within a subject).
# Neither carries an interaction, and the question a mixed design is usually run
# to answer -- does the within-subject effect DIFFER between the groups -- is
# exactly that interaction.
#
# Running the between-subjects kernel on a long file instead is worse than
# unavailable: a subject contributing two conditions appears as two rows, the
# permutation treats them as independent curves, and the test is
# anticonservative. The app already warns on import when identifiers repeat;
# this module is what that warning should point at.
#
# TWO KERNELS, and they answer different questions:
#
#   dance_mixed_fanova()   the whole curve. A smooth of time per cell of the
#                          design plus a per-subject random smooth, fitted by
#                          mgcv. The repeated measurement is IN the model rather
#                          than permuted around.
#
#   dance_mixed_cosinor()  the rhythm parameters. One linear mixed model over
#                          all observations, with the cosine/sine pair crossed
#                          with both factors and a random MESOR and random
#                          (cos, sin) per subject. This is a genuine mixed
#                          cosinor, NOT the app's two-stage route (fit per
#                          subject, then compare the estimates), which cannot
#                          respect pairing and treats every subject's estimate
#                          as if it were measured without error.
#
# Both are PURE: data frame in, list out. No `input`, no `values`, no
# notifications -- so they can be unit-tested, and emitted into the exported
# script by the deparse() machinery in server/90_export.R.
#
# WHAT THESE DO NOT DO. Neither is a permutation test. mgcv's and lme4's
# p-values are approximate: the smoothing parameters and the variance components
# were estimated from the same data, and the tests condition on those estimates.
# They are the standard tools for this design and their approximation is
# well-understood, but they are not the exact-permutation guarantee the one-way
# kernels give, and the readout says so rather than letting the reader assume.
# ==============================================================================

if (!exists("%||%", mode = "function")) `%||%` <- function(a, b) if (is.null(a)) b else a

# ------------------------------------------------------------------ long form
# Curves are held subjects x time. A mixed model needs one row per observation,
# and it needs the subject identifier to survive that reshaping -- which is the
# whole point, since that column is what makes the design mixed.
#
# `time_points` is the axis the curves were measured on. Pass real elapsed time
# when the columns are unevenly spaced: the smooth is a function of whatever is
# passed here, and a smooth of the column index on unevenly spaced data is a
# smooth of the wrong variable.
dance_mixed_long <- function(curves, time_points, subject, between, within) {
  curves <- as.matrix(curves)
  n <- nrow(curves); p <- ncol(curves)
  stopifnot(length(time_points) == p, length(subject) == n,
            length(between) == n, length(within) == n)

  d <- data.frame(
    .row     = rep(seq_len(n), each = p),
    subject  = factor(rep(as.character(subject), each = p)),
    between  = factor(rep(as.character(between), each = p)),
    within   = factor(rep(as.character(within),  each = p)),
    t        = rep(time_points, times = n),
    y        = as.vector(t(curves)),
    stringsAsFactors = FALSE
  )
  d <- d[is.finite(d$y), , drop = FALSE]
  # One smooth per cell of the design: mgcv's `by=` takes a single factor, so
  # the cell is formed here rather than by an interaction inside the formula.
  d$cell <- droplevels(interaction(d$within, d$between, sep = " x ", drop = TRUE))
  d$subject <- droplevels(d$subject)
  d$between <- droplevels(d$between)
  d$within  <- droplevels(d$within)
  d
}

# Is this actually a mixed design, and is it estimable? Returns a character
# vector of reasons it is not, empty when it is. Called before fitting so the
# refusal names the problem instead of surfacing as a fitting error.
dance_mixed_check <- function(d) {
  msg <- character(0)
  if (nlevels(d$within) < 2)
    msg <- c(msg, "The within-subject factor has fewer than 2 levels.")
  if (nlevels(d$between) < 2)
    msg <- c(msg, "The between-subject factor has fewer than 2 levels.")
  if (nlevels(d$subject) < 3)
    msg <- c(msg, "Fewer than 3 subjects: a random subject effect is not estimable.")

  # the between factor must not vary within a subject, or it is not between
  bt <- tapply(as.character(d$between), d$subject, function(x) length(unique(x)))
  if (any(bt > 1))
    msg <- c(msg, sprintf(
      "%d subject(s) have more than one level of the between-subject factor, so it is not a between-subject factor: %s.",
      sum(bt > 1), paste(names(bt)[bt > 1][seq_len(min(5, sum(bt > 1)))], collapse = ", ")))

  # the within factor must vary within a subject for at least some subjects
  wt <- tapply(as.character(d$within), d$subject, function(x) length(unique(x)))
  if (all(wt < 2))
    msg <- c(msg, "No subject has more than one level of the within-subject factor, so nothing is repeated within subject.")

  # every cell must be occupied, or the interaction is not identified
  tabc <- table(d$within, d$between)
  if (any(tabc == 0))
    msg <- c(msg, "At least one cell of the design is empty, so the interaction is not identified.")

  # AUDIT (P20/R8). The array builder writes Y[subject, condition, time] <- y by
  # index assignment, and index assignment keeps the LAST value written. So a
  # frame carrying two rows for the same (subject, condition, time) -- a
  # duplicated import, a repeated trial that was meant to be averaged, a file
  # concatenated with itself -- silently used whichever row happened to come
  # last, and REORDERING THE ROWS could then change every p-value in the result
  # while nothing anywhere said the data were ambiguous. The permutation schemes
  # assume one observation per design cell; that assumption is checked here
  # rather than assumed, and a violation names the offending keys instead of
  # reporting a count.
  key <- paste(as.character(d$subject), as.character(d$within), d$t, sep = "\r")
  dup <- unique(key[duplicated(key)])
  if (length(dup)) {
    show <- utils::head(dup, 5)
    parts <- vapply(strsplit(show, "\r", fixed = TRUE), function(p)
      sprintf("%s / %s / t = %s", p[1], p[2], p[3]), character(1))
    msg <- c(msg, sprintf(paste0(
      "%d participant-by-condition-by-time cell(s) appear more than once, so the ",
      "design has no single value per cell and the result would depend on the row ",
      "order of the file: %s%s. Aggregate the repeats (or drop the duplicate rows) ",
      "before running this."),
      length(dup), paste(parts, collapse = "; "),
      if (length(dup) > 5) sprintf(", and %d more", length(dup) - 5) else ""))
  }
  msg
}

# How balanced is it? Reported, never silently fixed: an unbalanced mixed design
# is analysable, but which subjects are missing which condition changes what the
# random effect can absorb, and the reader should be told.
dance_mixed_balance <- function(d) {
  wt <- tapply(as.character(d$within), d$subject, function(x) length(unique(x)))
  list(
    n_subjects   = nlevels(d$subject),
    n_complete   = sum(wt == nlevels(d$within)),
    n_partial    = sum(wt < nlevels(d$within)),
    cells        = as.data.frame(table(within = d$within, between = d$between),
                                 stringsAsFactors = FALSE),
    n_obs        = nrow(d),
    between_n    = vapply(split(d$subject, d$between),
                          function(s) length(unique(as.character(s))), integer(1))
  )
}

# ============================================================ mixed functional
# A smooth of time in each cell of the design, plus a factor-smooth random
# effect per subject. k_time bounds the per-cell smooth; k_subject the subject
# curves, which are deliberately coarser -- they are a nuisance term and giving
# them the same flexibility as the fixed smooths lets them absorb the effect.
#
# The interaction is tested by refitting without it and comparing: with
# penalised smooths there is no single clean F for "the cell smooths differ", so
# the comparison is what carries the claim.
dance_mixed_fanova <- function(d, k_time = 12, k_subject = 6, method = "fREML") {
  if (!requireNamespace("mgcv", quietly = TRUE))
    return(list(ok = FALSE, message = "mgcv is required for the mixed functional model."))
  bad <- dance_mixed_check(d)
  if (length(bad)) return(list(ok = FALSE, message = paste(bad, collapse = " ")))

  # k cannot exceed the number of distinct time points the data actually has
  n_t <- length(unique(d$t))
  k_time    <- max(3L, min(as.integer(k_time), n_t - 1L))
  k_subject <- max(3L, min(as.integer(k_subject), n_t - 1L))

  fit_one <- function(form) tryCatch(
    mgcv::bam(form, data = d, method = method, discrete = TRUE),
    error = function(e) tryCatch(mgcv::bam(form, data = d, method = method),
                                 error = function(e2) NULL))

  f_full <- stats::as.formula(sprintf(
    "y ~ within * between + s(t, by = cell, k = %d) + s(t, subject, bs = 'fs', k = %d, m = 1)",
    k_time, k_subject))
  f_add  <- stats::as.formula(sprintf(
    "y ~ within + between + s(t, by = within, k = %d) + s(t, by = between, k = %d) + s(t, subject, bs = 'fs', k = %d, m = 1)",
    k_time, k_time, k_subject))

  m_full <- fit_one(f_full)
  if (is.null(m_full)) return(list(ok = FALSE, message = "The mixed functional model did not converge."))
  m_add  <- fit_one(f_add)

  an <- tryCatch(mgcv::anova.gam(m_full), error = function(e) NULL)

  # AUDIT (P18.5). The comparison used to read AIC off the fREML fits. The two
  # models differ in their FIXED-effect structure -- the full one has the
  # within x between interaction and cell-specific smooths, the additive one
  # does not -- and a restricted likelihood is computed on a different set of
  # contrasts for each, so the two scores are not on a common scale. The
  # reported model keeps its fREML fit, because that is the better basis for
  # the smoothing parameters and hence for the fitted curves; the COMPARISON is
  # made on a pair refitted by ML.
  #
  # On the data this was raised against the verdict does not change (+33.5 by
  # fREML, +32.1 by ML), which is worth stating rather than hiding: the fix is
  # for the method, not because the answer was wrong.
  fit_ml <- function(form) tryCatch(
    mgcv::bam(form, data = d, method = "ML"), error = function(e) NULL)
  m_full_ml <- fit_ml(f_full)
  m_add_ml  <- fit_ml(f_add)
  aic_basis <- if (!is.null(m_full_ml) && !is.null(m_add_ml)) "ML" else method

  aic_full <- tryCatch(stats::AIC(m_full_ml %||% m_full), error = function(e) NA_real_)
  aic_add  <- {
    mm <- m_add_ml %||% m_add
    if (is.null(mm)) NA_real_ else tryCatch(stats::AIC(mm), error = function(e) NA_real_)
  }

  list(
    ok = TRUE,
    model        = m_full,
    model_add    = m_add,
    formula_full = deparse1(f_full),
    formula_add  = if (is.null(m_add)) NA_character_ else deparse1(f_add),
    s_table      = if (is.null(an)) NULL else as.data.frame(an$s.table),
    p_table      = if (is.null(an)) NULL else as.data.frame(an$p.table),
    aic_full     = aic_full,
    aic_additive = aic_add,
    aic_delta    = aic_add - aic_full,     # > 0 favours keeping the interaction
    aic_basis    = aic_basis,              # "ML" when the comparison refit worked
    dev_expl     = tryCatch(summary(m_full)$dev.expl, error = function(e) NA_real_),
    k_time       = k_time,
    k_subject    = k_subject,
    method       = method,
    balance      = dance_mixed_balance(d),
    n_obs        = nrow(d)
  )
}

# Fitted cell means over a time grid, for plotting and for the report. The
# subject term is excluded, so these are the POPULATION curves, not any
# particular subject's.
dance_mixed_fanova_curves <- function(res, n_grid = 100) {
  if (!isTRUE(res$ok) || is.null(res$model)) return(NULL)
  m <- res$model
  d <- m$model
  cells <- levels(d$cell)
  tr <- range(d$t, na.rm = TRUE)
  tg <- seq(tr[1], tr[2], length.out = n_grid)
  sub1 <- levels(d$subject)[1]
  nd <- do.call(rbind, lapply(cells, function(cl) {
    parts <- strsplit(cl, " x ", fixed = TRUE)[[1]]
    data.frame(t = tg, cell = factor(cl, levels = cells),
               within  = factor(parts[1], levels = levels(d$within)),
               between = factor(parts[2], levels = levels(d$between)),
               subject = factor(sub1, levels = levels(d$subject)),
               stringsAsFactors = FALSE)
  }))
  pr <- tryCatch(
    stats::predict(m, newdata = nd, se.fit = TRUE,
                   exclude = grep("subject", sapply(m$smooth, function(s) s$label),
                                  value = TRUE)),
    error = function(e) NULL)
  if (is.null(pr)) return(NULL)
  nd$fit <- as.numeric(pr$fit); nd$se <- as.numeric(pr$se.fit)
  nd
}

# ================================================================ mixed cosinor
# y ~ (cos + sin) * within * between + (1 + cos + sin | subject)
#
# The cosinor is linear in the (cos, sin) pair once the period is fixed, so a
# linear mixed model fits it directly: the fixed part gives a MESOR, amplitude
# and acrophase per cell of the design, and the random part lets each subject
# have their own MESOR and their own rhythm.
#
# Amplitude and acrophase are recovered from the cell's (cos, sin) coefficients.
# Their standard errors are NOT propagated here: amplitude is a norm and
# acrophase an angle, both nonlinear in the coefficients, and a delta-method SE
# on a phase near the boundary is misleading. What is tested instead is what the
# model can test exactly -- whether the (cos, sin) pair differs between cells,
# which is the rhythm differing, on 2 degrees of freedom.
dance_mixed_cosinor <- function(d, period = 24, n_harmonics = 1) {
  if (!requireNamespace("lme4", quietly = TRUE))
    return(list(ok = FALSE, message = "lme4 is required for the mixed cosinor."))
  bad <- dance_mixed_check(d)
  if (length(bad)) return(list(ok = FALSE, message = paste(bad, collapse = " ")))
  n_harmonics <- max(1L, as.integer(n_harmonics))

  dd <- d
  harm <- character(0)
  for (h in seq_len(n_harmonics)) {
    w <- 2 * pi * h / period
    dd[[paste0("c", h)]] <- cos(w * dd$t)
    dd[[paste0("s", h)]] <- sin(w * dd$t)
    harm <- c(harm, paste0("c", h), paste0("s", h))
  }
  hterms <- paste0("(", paste(harm, collapse = " + "), ")")
  f <- stats::as.formula(sprintf("y ~ %s * within * between + (1 + %s | subject)",
                                 hterms, paste(harm, collapse = " + ")))

  m <- tryCatch(
    lme4::lmer(f, data = dd, REML = TRUE,
               control = lme4::lmerControl(optimizer = "bobyqa",
                                           optCtrl = list(maxfun = 2e5))),
    error = function(e) NULL)
  # A full random rhythm can be too much for a small sample. Fall back to a
  # random MESOR only, and SAY which was fitted -- the two are different models
  # and a reader cannot tell them apart from the coefficients.
  simplified <- FALSE
  if (is.null(m)) {
    f2 <- stats::as.formula(sprintf("y ~ %s * within * between + (1 | subject)", hterms))
    m <- tryCatch(lme4::lmer(f2, data = dd, REML = TRUE,
                             control = lme4::lmerControl(optimizer = "bobyqa")),
                  error = function(e) NULL)
    simplified <- !is.null(m)
    f <- f2
  }
  if (is.null(m)) return(list(ok = FALSE, message = "The mixed cosinor did not converge."))

  fe <- lme4::fixef(m)
  cells <- expand.grid(within = levels(dd$within), between = levels(dd$between),
                       stringsAsFactors = FALSE)
  ref_w <- levels(dd$within)[1]; ref_b <- levels(dd$between)[1]
  getc <- function(nm) if (nm %in% names(fe)) unname(fe[[nm]]) else 0

  # Sum the terms that are active in a cell. Written out rather than taken from
  # predict(), because the amplitude and acrophase need the (cos, sin) PAIR for
  # that cell, not a fitted value.
  coef_in_cell <- function(base, w, b) {
    v <- getc(base)
    if (w != ref_w) v <- v + getc(paste0(base, ":within", w))
    if (b != ref_b) v <- v + getc(paste0(base, ":between", b))
    if (w != ref_w && b != ref_b)
      v <- v + getc(paste0(base, ":within", w, ":between", b))
    v
  }
  mesor_in_cell <- function(w, b) {
    v <- getc("(Intercept)")
    if (w != ref_w) v <- v + getc(paste0("within", w))
    if (b != ref_b) v <- v + getc(paste0("between", b))
    if (w != ref_w && b != ref_b) v <- v + getc(paste0("within", w, ":between", b))
    v
  }

  rows <- do.call(rbind, lapply(seq_len(nrow(cells)), function(i) {
    w <- cells$within[i]; b <- cells$between[i]
    out <- data.frame(within = w, between = b, mesor = mesor_in_cell(w, b),
                      stringsAsFactors = FALSE)
    for (h in seq_len(n_harmonics)) {
      bc <- coef_in_cell(paste0("c", h), w, b)
      bs <- coef_in_cell(paste0("s", h), w, b)
      amp <- sqrt(bc^2 + bs^2)
      # acrophase in TIME units, on the effective period of harmonic h
      acro <- ((atan2(bs, bc)) %% (2 * pi)) / (2 * pi * h / period)
      out[[paste0("amplitude_", h)]] <- amp
      out[[paste0("acrophase_", h)]] <- acro
    }
    out
  }))

  # Does the rhythm differ? A likelihood-ratio test on the (cos, sin) pair,
  # which is the 2-df question "is there any difference in amplitude OR phase",
  # refitted by ML because REML likelihoods are not comparable across fixed
  # effects.
  # Terms are dropped with update(), not rebuilt with reformulate(): the random
  # part appears in term.labels as the string "1 + c1 + s1 | subject", and
  # reformulate() parses text, so putting that back through it produces a
  # different model. update(m, . ~ . - <terms>) leaves the random part alone.
  #
  # Marginality is respected: a test of the within factor on the rhythm drops
  # its higher-order terms too, so it asks "does the within factor affect the
  # rhythm AT ALL", not "does it affect it after allowing it to affect it".
  hx <- function(suffix) paste0(harm, ":", suffix)     # c1:within, s1:within, ...
  m_ml <- tryCatch(stats::update(m, REML = FALSE), error = function(e) NULL)

  lrt <- function(drop_terms, label) {
    if (is.null(m_ml)) return(NULL)
    drop_terms <- intersect(drop_terms, attr(stats::terms(f), "term.labels"))
    if (!length(drop_terms)) return(NULL)
    delta <- stats::as.formula(paste(". ~ . -", paste(drop_terms, collapse = " - ")))
    m0 <- tryCatch(stats::update(m_ml, delta), error = function(e) NULL)
    if (is.null(m0)) return(NULL)
    a <- tryCatch(stats::anova(m0, m_ml), error = function(e) NULL)
    if (is.null(a) || nrow(a) < 2) return(NULL)
    data.frame(term = label, df = a$Df[2], chisq = a$Chisq[2],
               p = a$`Pr(>Chisq)`[2], stringsAsFactors = FALSE)
  }
  three <- hx("within:between")
  tests <- do.call(rbind, Filter(Negate(is.null), list(
    lrt(three, "rhythm: within x between interaction"),
    lrt(c(hx("within"), three), "rhythm differs by the within factor"),
    lrt(c(hx("between"), three), "rhythm differs by the between factor")
  )))

  list(
    ok = TRUE, model = m, formula = deparse1(f),
    period = period, n_harmonics = n_harmonics,
    cells = rows, tests = tests,
    random_rhythm = !simplified,
    n_obs = nrow(dd), balance = dance_mixed_balance(d),
    singular = tryCatch(lme4::isSingular(m), error = function(e) NA)
  )
}
# AUDIT (P19). One description of a fitted mixed model, used by every place that
# shows one -- the Functional ANOVA tab and the Cosinor tab both reach this
# module now, and two entry points writing their own readout is how the same fit
# comes to be described two different ways. This is the sixth time this audit has
# extracted a duplicated definition; the pattern is always the same.
dance_mixed_readout <- function(res) {

  if (is.null(res)) {
  cat("Run a mixed analysis to see results.\n\n")
  cat("This tab is for a design with one BETWEEN-subject factor and one\n")
  cat("WITHIN-subject factor. The Functional ANOVA tab handles one or the\n")
  cat("other, not both, and carries no interaction term.\n")
  return(invisible(NULL))
}
f2 <- function(x) if (is.finite(x)) sprintf("%.2f", x) else "--"
f3 <- function(x) if (is.finite(x)) sprintf("%.3f", x) else "--"
pf <- function(p) if (!is.finite(p)) "--" else if (p < .001) "< .001" else sub("^0", "", sprintf("%.3f", p))

b <- res$balance
cat("Design\n======\n")
cat(sprintf("  between: %s   within: %s   subjects: %d   observations: %d\n",
            res$between_name, res$within_name, b$n_subjects, res$n_obs))
cat(sprintf("  time axis: %s\n", res$time_axis))
if (b$n_partial > 0)
  cat(sprintf("  %d participant(s) do not have every level of the within factor.\n",
              b$n_partial))
cat("\n")

if (identical(res$kind, "permutation")) {
  # AUDIT (P20/R5). This heading used to read "exact permutation" and the note
  # below it asserted that ALL THREE effects were exact. That was wrong for two
  # of them and the error was not small: with 5 against 30 participants and a
  # 5:1 dispersion ratio in the within-subject contrast, the interaction test
  # rejected 75% of the time at a nominal 5%. Exchangeability requires the group
  # distributions to be IDENTICAL, not merely to have equal means, and unequal
  # group sizes with unequal dispersion is the classic Behrens-Fisher situation
  # in which relabelling is not a symmetry. The statistic is studentised now,
  # which fixes most of it, and each effect states its own status rather than
  # borrowing the within effect's.
  cat("Mixed functional ANOVA, permutation\n")
  cat("===================================\n")
  cat(sprintf("  %d participants x %d conditions; groups: %s\n",
              res$n_subjects, res$n_conditions,
              paste(sprintf("%s = %d", res$group_names, res$group_sizes), collapse = ", ")))
  cat(sprintf("  B = %d permutations; smallest attainable p = %.4f\n",
              res$n_permutations, res$p_floor))
  cat(sprintf("  pointwise p adjusted by %s at alpha = %s\n",
              res$correction, format(res$alpha)))
  if (!is.null(res$multiplicity_family))
    cat(sprintf("  multiplicity family: %s\n", res$multiplicity_family))
  cat("\n")
  lab <- c(within = "within-subject effect", between = "between-subject effect",
           interaction = "interaction")
  status_lab <- c(exact = "exact", calibrated = "asymptotic",
                  liberal = "global p ANTI-CONSERVATIVE")
  for (nm in c("within", "between", "interaction")) {
    r <- res[[nm]]
    if (is.null(r)) next
    st <- r$calibration$status %||% "unknown"
    cat(sprintf("  %-22s global p %s   significant at %d of %d points (%.1f%%)   [%s]\n",
                lab[[nm]], pf(r$global_p), r$n_significant,
                length(r$p_values), 100 * r$n_significant / length(r$p_values),
                status_lab[[st]] %||% st))
  }
  cat("\nThe permutation schemes\n-----------------------\n")
  cat("  within       condition labels relabelled INSIDE each participant\n")
  cat("  between      whole participants relabelled across groups\n")
  cat("  interaction  group labels relabelled on the subject-centred profiles\n")
  cat("\nWhat each one establishes\n-------------------------\n")
  for (nm in c("within", "between", "interaction")) {
    r <- res[[nm]]
    if (is.null(r) || is.null(r$calibration)) next
    cat(sprintf("  %s:\n", lab[[nm]]))
    cat(strwrap(r$calibration$message, width = 72, prefix = "    "), sep = "\n")
    cat("\n")
  }
  cat("The interaction scheme itself is the one usually said not to exist: an\n")
  cat("interaction is a between-group difference in the within-subject contrast,\n")
  cat("and under that null the subject-centred profiles are exchangeable across\n")
  cat("groups whatever the two main effects do. What the studentisation adds is\n")
  cat("validity when those profiles have DIFFERENT DISPERSIONS in the two groups,\n")
  cat("which exchangeability alone does not cover. Calibration is measured in\n")
  cat("tests/mixed_calibration_test.R.\n")
  cat("\nWhat this does not establish\n----------------------------\n")
  cat("  A pointwise procedure answers WHERE the curves differ, not whether\n")
  cat("  they differ overall; the global statistic is the whether. The\n")
  cat("  correction controls the expected proportion of false positives among\n")
  cat("  the flagged points, not the familywise error. And no permutation p\n")
  cat("  can fall below 1/(B+1), so identical p values in a table are the\n")
  cat("  resolution limit, not a tie.\n")

} else if (identical(res$kind, "fanova") || identical(res$kind, "model")) {
  cat("Mixed functional model\n======================\n")
  cat("  ", res$formula_full, "\n\n", sep = "")
  if (!is.null(res$s_table)) {
    cat("Smooth terms (approximate)\n")
    st <- res$s_table
    for (i in seq_len(nrow(st)))
      cat(sprintf("  %-28s edf %7s   F %8s   p %s\n", rownames(st)[i],
                  f2(st[i, "edf"]), f2(st[i, "F"]), pf(st[i, "p-value"])))
    cat("\n")
  }
  if (!is.null(res$p_table)) {
    cat("Parametric terms (mean level)\n")
    pt <- res$p_table
    for (i in seq_len(nrow(pt)))
      cat(sprintf("  %-28s est %8s   SE %7s   t %7s   p %s\n", rownames(pt)[i],
                  f2(pt[i, "Estimate"]), f2(pt[i, "Std. Error"]),
                  f2(pt[i, "t value"]), pf(pt[i, "Pr(>|t|)"])))
    cat("\n")
  }
  cat(sprintf("Deviance explained: %s%%\n", f2(100 * res$dev_expl)))
  if (is.finite(res$aic_delta)) {
    cat(sprintf("Interaction in the SHAPE (%s-based model comparison)\n",
                res$aic_basis %||% "ML"))
    cat(sprintf("  AIC %s with the interaction, %s without; delta %+.1f\n",
                f2(res$aic_full), f2(res$aic_additive), res$aic_delta))
    # P18.5: reported as model-comparison evidence, not as a hypothesis test.
    # A delta is a weight of evidence; calling it a decision at 2 units dresses
    # a continuous quantity as a verdict.
    cat(sprintf("  %s\n", if (res$aic_delta > 2)
      "the comparison favours letting each cell have its own temporal shape"
      else if (res$aic_delta < -2)
      "the comparison favours the additive model: no support for cell-specific shapes"
      else "the two models are within 2 AIC: this comparison does not separate them"))
    cat("  This is evidence for one model over another, not a test of a null\n")
    cat("  hypothesis, and no p-value should be quoted from it.\n")
  }
  cat("\nWhat this does not establish\n----------------------------\n")
  cat("  The p-values are APPROXIMATE. The smoothing parameters were estimated\n")
  cat("  from these data and the tests condition on those estimates; this is not\n")
  cat("  the exact permutation guarantee the one-way fANOVA gives. The subject\n")
  cat("  term is a random functional effect, so the cell smooths are population\n")
  cat("  curves and no single participant's curve is claimed. A smooth term's\n")
  cat("  p-value tests whether that cell's curve is flat, NOT whether two cells\n")
  cat("  differ -- the interaction comparison above is what addresses that.\n")

} else {
  cat("Mixed cosinor\n=============\n")
  cat("  ", res$formula, "\n", sep = "")
  cat(sprintf("  period %s, %d harmonic(s); random rhythm per subject: %s\n\n",
              format(res$period), res$n_harmonics,
              if (isTRUE(res$random_rhythm)) "yes" else
                "NO -- the full model did not converge, so only a random MESOR was fitted"))
  cl <- res$cells
  cat(sprintf("  %-14s %-10s %8s %10s %10s\n", res$within_name, res$between_name,
              "MESOR", "amplitude", "acrophase"))
  for (i in seq_len(nrow(cl)))
    cat(sprintf("  %-14s %-10s %8s %10s %10s\n", cl$within[i], cl$between[i],
                f2(cl$mesor[i]), f2(cl$amplitude_1[i]), f2(cl$acrophase_1[i])))
  cat("\n  Acrophase is in time units from the start of the observation window.\n")
  if (!is.null(res$tests)) {
    cat("\nLikelihood-ratio tests on the (cosine, sine) pair\n")
    for (i in seq_len(nrow(res$tests)))
      cat(sprintf("  %-40s chi2(%d) = %s, p %s\n", res$tests$term[i],
                  res$tests$df[i], f2(res$tests$chisq[i]), pf(res$tests$p[i])))
    cat("\n  Each test drops a cosine/sine pair, so it asks whether the rhythm\n")
    cat("  differs in amplitude OR phase -- a 2-df question per harmonic. A test\n")
    cat("  of one factor also drops its higher-order terms, so it is the effect\n")
    cat("  of that factor overall, not conditional on the interaction.\n")
  }
  if (isTRUE(res$singular))
    cat("\n  WARNING: the random-effects fit is SINGULAR. A variance component is\n",
        "  estimated at zero; the rhythm parameters are still usable but the\n",
        "  random structure is not supported by this sample.\n", sep = "")
  cat("\nWhat this does not establish\n----------------------------\n")
  cat("  No standard errors are quoted for amplitude or acrophase. Amplitude is\n")
  cat("  a norm and acrophase an angle, both nonlinear in the coefficients, and\n")
  cat("  a delta-method interval on a phase near the period boundary misleads.\n")
  cat("  The tests above are on the coefficient pair, which is what the model\n")
  cat("  can test exactly. The period was FIXED, not estimated, so everything\n")
  cat("  here is conditional on that choice.\n")
}

  invisible(NULL)
}

