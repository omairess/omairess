# ==============================================================================
# tests/traj_framework_test.R — P21 phase 2: the trajectory framework, layers A-D
# ==============================================================================
# The framework replaces a two-stage procedure whose group test was a t-test on
# per-participant point estimates (finding A1) and whose mixed path was
# hard-coded to exactly one between x one within factor (finding A6). What has
# to be shown is therefore not only that it runs, but that it is GENERAL -- the
# same call takes 2, 3 or 4 factors of either kind -- and that it recovers what
# was put in.
#
# Simulation counts here are deliberately modest and the bands around them are
# correspondingly wide: this is a correctness and generality suite, not the
# full §20 calibration grid, which belongs to the validation gate between
# phases 2 and 3. DANCE_TRAJ_NSIM raises the counts.
#
# Run with:  Rscript tests/traj_framework_test.R      (from the DANCE dir)
# ==============================================================================
.libPaths(c("~/Rlib", .libPaths()))
suppressMessages({ library(lme4); library(lmerTest); library(emmeans) })
`%||%` <- function(a, b) if (is.null(a)) b else a
ok_n <- 0L; bad <- 0L
chk <- function(cond, good, bad_msg) {
  if (isTRUE(cond)) { cat("ok   ", good, "\n"); ok_n <<- ok_n + 1L }
  else { cat("FAIL:", bad_msg, "\n"); bad <<- bad + 1L }
}
app_dir <- if (dir.exists("server")) "." else "DANCE"
e <- new.env(parent = globalenv())
for (f in c("server/08d_helpers_traj.R", "server/08e_helpers_trajfit.R",
            "server/08f_helpers_trajinf.R"))
  source(file.path(app_dir, f), local = e)

NSIM <- suppressWarnings(as.integer(Sys.getenv("DANCE_TRAJ_NSIM")))
if (is.na(NSIM) || NSIM < 1) NSIM <- 40L

# ------------------------------------------------------------------ generator
# `shift` is applied in RADIANS to the named cell only, so the truth is known.
make <- function(seed, n_per_group = 6, groups = c("a", "b"),
                 conds = c("p", "q"), nt = 12, amp = 4, shift = 0,
                 shift_cell = NULL, level = 0, level_group = NULL,
                 noise = 1, slope = 0) {
  set.seed(seed)
  tp <- seq(0, 22, length.out = nt)
  subj <- sprintf("S%03d", seq_len(n_per_group * length(groups)))
  sg <- rep(groups, each = n_per_group)
  rows <- list(); k <- 1L
  meta <- list(subject = character(0), Group = character(0), Condition = character(0))
  # TWO LEVELS OF RANDOM VARIATION, and the first draft had only the second.
  #
  # A participant offset (sb, sa, sp) is shared by both of that participant's
  # curves; a curve offset is drawn afresh for each condition. That is what
  # repeated-measures data looks like: the same person is rhythmic in a way that
  # persists across conditions AND differs between them. The first draft drew
  # base/a/ph once per CURVE and nothing at participant level, which makes two
  # curves from one person no more alike than two from different people -- a
  # null the app will never meet, and one that happens to be easier on the model
  # than the real thing. Adding the participant level makes this generator
  # strictly harder, not kinder: it is the case whose type-I error measured
  # 0.283 before the curve-level rung existed.
  for (i in seq_along(subj)) {
    sb <- rnorm(1, 0, 2); sa <- rnorm(1, 0, .4); sp <- rnorm(1, 0, .2)
    for (cc in conds) {
      base <- 10 + sb + rnorm(1, 0, 2)
      a <- amp + sa + rnorm(1, 0, .4)
      ph <- 2 + sp + rnorm(1, 0, .2)
      if (!is.null(shift_cell) && paste(sg[i], cc) == shift_cell) ph <- ph + shift
      lv <- if (!is.null(level_group) && sg[i] == level_group) level else 0
      rows[[k]] <- base + lv + a * cos(2 * pi * tp / 24 - ph) + slope * tp +
                   rnorm(nt, 0, noise)
      meta$subject <- c(meta$subject, subj[i])
      meta$Group <- c(meta$Group, sg[i])
      meta$Condition <- c(meta$Condition, cc)
      k <- k + 1L
    }
  }
  list(Y = do.call(rbind, rows), t = tp, meta = meta)
}
build <- function(g, ...) e$dance_traj_long(g$Y, g$t, g$meta$subject,
                                            list(Group = g$meta$Group,
                                                 Condition = g$meta$Condition))

cat("-- layer A: the design is read off the data, not declared ------------\n")
g <- make(1); d <- build(g)
cls <- e$dance_traj_classify(d)
chk(cls$role[cls$factor == "Group"] == "between" &&
      cls$role[cls$factor == "Condition"] == "within",
    "Group is classified between-participant and Condition within-participant",
    sprintf("classification wrong: %s", paste(cls$role, collapse = "/")))
chk(identical(e$dance_traj_design_kind(cls), "mixed (between x within)"),
    "the design is named mixed (between x within)",
    "the design kind is wrong")
# a factor that is within for some participants and between for others
d2 <- d; d2$Condition[d2$subject == levels(d2$subject)[1]] <- "p"
cls2 <- e$dance_traj_classify(d2)
chk(cls2$role[cls2$factor == "Condition"] == "partial",
    "a factor that is within for some participants and between for others is 'partial'",
    "an inconsistent factor was classified as clean")
sp_bad <- e$dance_traj_spec(d2)
chk(isFALSE(sp_bad$ok) && grepl("neither a between", sp_bad$message),
    "and the spec refuses it rather than guessing",
    "a partial factor was accepted into the model")
# a numeric covariate is not a design factor
d3 <- e$dance_traj_long(g$Y, g$t, g$meta$subject,
                        list(Group = g$meta$Group, Condition = g$meta$Condition,
                             Age = round(runif(nrow(g$Y), 20, 70))))
cls3 <- e$dance_traj_classify(d3)
chk(cls3$role[cls3$factor == "Age"] == "covariate",
    "a numeric covariate is classified as a covariate, not a factor",
    "a numeric covariate was treated as a design factor")

cat("\n-- layer A: the formula is not specific to any design size ----------\n")
f2 <- e$dance_traj_formula(c("c1", "s1"), c("Group", "Condition"))
chk(identical(f2, "y ~ (c1 + s1) * Group * Condition"),
    "2 factors: y ~ (c1 + s1) * Group * Condition",
    sprintf("got %s", f2))
f3 <- e$dance_traj_formula(c("trend_lin", "c1", "s1", "c2", "s2"),
                           c("Group", "Condition", "Site"))
chk(identical(f3, "y ~ (trend_lin + c1 + s1 + c2 + s2) * Group * Condition * Site"),
    "3 factors, 2 harmonics and a trend: the same call, a longer term list",
    sprintf("got %s", f3))
f4 <- e$dance_traj_formula(c("c1", "s1"), c("Group"), covariates = "Age")
chk(grepl("+ Age", f4, fixed = TRUE) && !grepl("Age *", f4, fixed = TRUE),
    "a covariate enters additively, not crossed with the whole harmonic block",
    sprintf("covariate handling wrong: %s", f4))
f0 <- e$dance_traj_formula(c("c1", "s1"), character(0))
chk(identical(f0, "y ~ (c1 + s1)"),
    "no design factors: one trajectory for the whole sample",
    sprintf("got %s", f0))

cat("\n-- layer A: basis, reference time and conditional linearity ---------\n")
b <- e$dance_traj_basis(d, period = 24, n_harmonics = 2, trend = "linear")
chk(identical(b$harm_terms, c("c1", "s1", "c2", "s2")) &&
      identical(b$trend_terms, "trend_lin"),
    "two harmonics give c1, s1, c2, s2 alongside the trend column",
    "the basis terms are wrong")
chk(abs(b$t0 - min(d$t)) < 1e-12 && abs(b$data$trend_lin[1]) < 1e-12,
    sprintf("t = 0 is set at the first observation (%.3f) and the trend column starts there",
            b$t0),
    "the reference time is not where it claims to be")
chk(max(abs(b$data$c1 - cos(2 * pi * (b$data$t - b$t0) / 24))) < 1e-12,
    "the first harmonic is cos(2*pi*(t - t0)/P) exactly",
    "the harmonic column is not what it says")
# exp_sat is linear once tau is fixed -- the property the tau profile rests on
bs <- e$dance_traj_basis(d, 24, 1, "exp_sat", tau = 8)
chk(max(abs(bs$data$trend_sat - (1 - exp(-(bs$data$t - bs$t0) / 8)))) < 1e-12,
    "at fixed tau the saturating trend is an ordinary fixed basis column",
    "the saturating column is wrong")
chk(inherits(try(e$dance_traj_basis(d, 24, 1, "exp_sat"), silent = TRUE), "try-error"),
    "and it refuses to build without a tau, rather than inventing one",
    "exp_sat accepted a missing tau")

cat("\n-- layer A: the random-effects ladder -------------------------------\n")
lad <- e$dance_traj_re_ladder(c("trend_lin", "c1", "s1"), c("c1", "s1"), "Condition")
chk(length(lad) == 7 && identical(lad[[length(lad)]]$formula, "(1 | subject)"),
    sprintf("%d rungs, most complex first, ending at (1 | subject)", length(lad)),
    sprintf("the ladder is wrong: %d rungs", length(lad)))
# The top rung must give each CURVE its own rhythm, not merely each participant
# an offset between conditions -- that distinction is the whole type-I fix.
chk(grepl("(1 + trend_lin + c1 + s1 | subject:Condition)", lad[[1]]$formula, fixed = TRUE) &&
      grepl("(1 + trend_lin + c1 + s1 | subject)", lad[[1]]$formula, fixed = TRUE),
    "the ladder starts at a curve-specific rhythm, not a within-factor offset",
    sprintf("the top rung is %s", lad[[1]]$formula))
lad_nw <- e$dance_traj_re_ladder(c("c1", "s1"), c("c1", "s1"))
chk(!any(grepl(":", vapply(lad_nw, function(x) x$formula, character(1)), fixed = TRUE)),
    "with no within-participant factor there is no curve level, and no such rung",
    "a curve-level rung appeared in a purely between-participant design")
lad0 <- e$dance_traj_re_ladder(c("c1", "s1"), c("c1", "s1"))
chk(!any(duplicated(vapply(lad0, function(x) x$formula, character(1)))),
    "structurally identical rungs are de-duplicated when there is no trend",
    "the ladder repeats a rung")

cat("\n-- layer A: identifiability refusals --------------------------------\n")
# 4 distinct times cannot support 3 harmonics plus a trend (8 basis columns),
# however many participants contribute -- the rank is set by the time grid.
g_small <- make(9, n_per_group = 3, nt = 4)
d_small <- build(g_small)
sp_thin <- e$dance_traj_spec(d_small, 24, 3, "linear")
chk(isFALSE(sp_thin$ok) && grepl("distinct time points", sp_thin$message),
    "3 harmonics on a 4-point time grid is refused for rank, not for row count",
    sprintf("a rank-deficient model was accepted: %s", sp_thin$message %||% "(no message)"))
chk(isFALSE(sp_thin$ok) && grepl("at most 1 harmonic", sp_thin$message),
    "and the refusal says how many harmonics the grid WOULD support",
    "the refusal does not say what would fit")
# but the same grid takes one harmonic with no trend
chk(isTRUE(e$dance_traj_spec(d_small, 24, 1, "none")$ok),
    "one harmonic with no trend on the same grid is accepted",
    "a model that fits within the rank was refused")

cat("\n-- layer B: the ladder is walked and the rung REPORTED --------------\n")
sp <- e$dance_traj_spec(d, 24, 1, "linear")
stopifnot(isTRUE(sp$ok))
fit <- e$dance_traj_fit(sp)
chk(isTRUE(fit$ok), sprintf("the model fits, at rung %d of %d (%s)",
                            fit$re_rung, fit$n_rungs, fit$re_label),
    "the model did not fit at any rung")
chk(identical(fit$simplified, fit$re_rung > 1L) &&
      (!fit$simplified || grepl("DIFFERENT model", fit$note)),
    "a simplification is reported as a different model, never silently",
    "a simplified fit did not say so")
chk(length(fit$attempts) == fit$re_rung &&
      all(vapply(fit$attempts[-fit$re_rung], function(a) !a$fitted || !a$converged, logical(1))),
    "every rung that was tried and rejected is logged, and only non-convergence rejects one",
    "the attempt log does not account for the rungs tried")
# Singularity must NOT descend the ladder: on the null data used to calibrate
# layer C the top rung was singular about half the time, and dropping it there
# is exactly what inflated the type-I error to 0.200.
# Data with NO participant-level variance at all: the participant terms sit on
# the boundary by construction, so the top rung is singular and must be kept.
gs <- make(77, n_per_group = 6)
gs$Y <- gs$Y - rowMeans(gs$Y) + 10        # every curve the same mean, no level spread
fit_s <- e$dance_traj_fit(e$dance_traj_spec(build(gs), 24, 1, "none"))
chk(isTRUE(fit_s$ok) && isTRUE(fit_s$singular) && fit_s$re_rung == 1L,
    "a singular top rung is KEPT, not descended past",
    sprintf("singular handling wrong: ok=%s rung=%s singular=%s",
            fit_s$ok, fit_s$re_rung %||% NA, fit_s$singular %||% NA))
chk(grepl("boundary", fit_s$note %||% "", fixed = TRUE) &&
      grepl("KEPT", fit_s$note %||% "", fixed = TRUE),
    "and the boundary is reported in the note rather than left for the reader to infer",
    "a singular fit came back without saying so")
chk(isFALSE(e$dance_traj_fit(sp, engine = "lmer", ar1 = TRUE)$ok),
    "a residual AR(1) request on lme4 is refused, naming glmmTMB",
    "lme4 silently accepted a residual correlation structure")

cat("\n-- layer C: the omnibus is a BLOCK test -----------------------------\n")
o <- e$dance_traj_omnibus(fit, "circadian")
chk(isTRUE(o$ok) && length(o$terms) == 6 &&
      all(c("c1:Group", "s1:Group", "c1:Group:Condition") %in% o$terms),
    sprintf("the circadian block is all 6 harmonic x design terms, tested jointly"),
    "the circadian block is not the full harmonic x design set")
chk(grepl("Kenward-Roger", o$method),
    sprintf("%d participants (24 curves) -> %s", sp$n_participants, o$method),
    sprintf("the wrong df method was chosen: %s", o$method))
ot <- e$dance_traj_omnibus(fit, "trajectory")
chk(isTRUE(ot$ok) && all(o$terms %in% ot$terms) && length(ot$terms) > length(o$terms),
    "the trajectory block strictly contains the circadian block",
    "the two blocks are not nested as they must be")
ol <- e$dance_traj_omnibus(fit, "level")
chk(isTRUE(ol$ok) && !any(grepl("c1|s1|trend", ol$terms)),
    "the level block holds no basis term -- it is the value at t = 0, not a MESOR",
    "the level block is contaminated with basis terms")

cat("\n-- layer C: it finds what was planted, and not what was not ---------\n")
g_ph <- make(21, n_per_group = 8, shift = 1.2, shift_cell = "b q")
f_ph <- e$dance_traj_fit(e$dance_traj_spec(build(g_ph), 24, 1, "none"))
o_ph <- e$dance_traj_omnibus(f_ph, "circadian")
chk(isTRUE(o_ph$ok) && o_ph$p < .01,
    sprintf("a planted 1.2 rad phase shift in one cell is detected (p = %.3g)", o_ph$p),
    sprintf("a planted phase shift was missed (p = %s)", format(o_ph$p)))
g_lv <- make(22, n_per_group = 8, level = 6, level_group = "b")
f_lv <- e$dance_traj_fit(e$dance_traj_spec(build(g_lv), 24, 1, "none"))
chk(e$dance_traj_omnibus(f_lv, "level")$p < .01 &&
      e$dance_traj_omnibus(f_lv, "circadian")$p > .05,
    "a pure level difference shows in the level block and NOT in the circadian one",
    "a level difference leaked into the circadian block")

cat(sprintf("\n-- layer C: size AND power under a true null (%d simulations each) ---\n", NSIM))
# THIS CELL ONCE FAILED, AT 0.275 AGAINST A NOMINAL .05, AND LAYER C SHIPPED
# UNVALIDATED BECAUSE OF IT. Two things were wrong with the random-effects
# ladder, both in server/08d and server/08e and both documented there:
#
#   1. the ladder had no curve (participant x condition) level at all, so the
#      variance a repeated-measures design puts there went into the residual and
#      shrank the standard errors of the terms under test;
#   2. singularity descended the ladder anyway -- first through an explicit
#      isSingular test, then, once that was removed, through lme4 filing
#      "boundary (singular) fit" in the same message list it uses for genuine
#      convergence failures. Between them the curve-level rung was dropped in
#      64-100% of null fits, i.e. almost always.
#
# SIZE IS ASSERTED WITH POWER BESIDE IT, on purpose. A size assertion alone is
# satisfied by a test that never rejects anything, and the fix here deliberately
# keeps a random structure whose components often sit on the boundary -- exactly
# the change that could buy calibration with power. Both numbers are measured on
# every run, so neither can rot into a comment.
null_p <- vapply(seq_len(NSIM), function(i) {
  ff <- e$dance_traj_fit(e$dance_traj_spec(build(make(3000 + i, n_per_group = 8)), 24, 1, "none"))
  if (!isTRUE(ff$ok)) return(NA_real_)
  e$dance_traj_omnibus(ff, "circadian")$p %||% NA_real_
}, numeric(1))
null_p <- null_p[is.finite(null_p)]
rate <- mean(null_p <= .05)
band <- 3 * sqrt(.05 * .95 / length(null_p))
chk(rate <= .05 + band,
    sprintf("the circadian block holds its size: %.3f at a nominal .05 (n = %d, 3-SE band <= %.3f)",
            rate, length(null_p), .05 + band),
    sprintf("the circadian block rejects %.3f of true nulls at a nominal .05 (n = %d), above the %.3f band",
            rate, length(null_p), .05 + band))

pow_p <- vapply(seq_len(NSIM), function(i) {
  g <- make(5000 + i, n_per_group = 8, shift = 0.5, shift_cell = "a q")
  ff <- e$dance_traj_fit(e$dance_traj_spec(build(g), 24, 1, "none"))
  if (!isTRUE(ff$ok)) return(NA_real_)
  e$dance_traj_omnibus(ff, "circadian")$p %||% NA_real_
}, numeric(1))
pow_p <- pow_p[is.finite(pow_p)]
pow <- mean(pow_p <= .05)
chk(pow >= .60,
    sprintf("and it was not bought with power: %.3f against a 0.5 rad phase shift in one cell (n = %d)",
            pow, length(pow_p)),
    sprintf("power against a 0.5 rad shift has fallen to %.3f (n = %d) -- the size result may be a dead test",
            pow, length(pow_p)))

# The result must still SAY which grid it was calibrated on, and must refuse to
# borrow that warrant outside it.
ff_one <- e$dance_traj_fit(e$dance_traj_spec(build(make(3001, n_per_group = 8)), 24, 1, "none"))
o_one <- e$dance_traj_omnibus(ff_one, "circadian")
chk(isTRUE(o_one$validated) && grepl("CALIBRATED", o_one$calibration %||% "", fixed = TRUE) &&
      grepl("0.040", o_one$calibration, fixed = TRUE),
    "a one-harmonic Kenward-Roger fit declares itself calibrated and names the measured rates",
    "the omnibus does not report its calibration grid")
ff_h2 <- e$dance_traj_fit(e$dance_traj_spec(build(make(3001, n_per_group = 8)), 24, 2, "none"))
o_h2 <- e$dance_traj_omnibus(ff_h2, "circadian")
chk(isFALSE(o_h2$validated) &&
      grepl("OUTSIDE THE CALIBRATED GRID", o_h2$calibration %||% "", fixed = TRUE) &&
      grepl("2 harmonics", o_h2$calibration, fixed = TRUE),
    "and a two-harmonic fit does NOT borrow it -- it says which way it left the grid",
    "a fit outside the simulated grid claimed the grid's calibration")

cat("\n-- layer D: parameters come from the SAME fit -----------------------\n")
co <- e$dance_traj_cell_coefs(fit, 1)
chk(isTRUE(co$ok) && length(co$cells) == 4 && abs(co$effective_period - 24) < 1e-12,
    "one (cos, sin) pair per design cell, on the harmonic's effective period",
    "the cell coefficients are wrong")
ap <- e$dance_traj_amp_phase(co)
chk(isTRUE(ap$ok) && all(abs(ap$table$amplitude - 4) < 1),
    sprintf("amplitudes recover the planted 4 (range %.2f-%.2f)",
            min(ap$table$amplitude), max(ap$table$amplitude)),
    "the amplitudes do not recover the truth")
chk(all(ap$table$amplitude_lo < ap$table$amplitude &
        ap$table$amplitude < ap$table$amplitude_hi),
    "every amplitude interval contains its own point estimate",
    "an amplitude interval does not contain its estimate")
# Bingham's rule: no rhythm at all -> the phase interval must be UNDEFINED
g_flat <- make(31, n_per_group = 8, amp = 0, noise = 3)
f_flat <- e$dance_traj_fit(e$dance_traj_spec(build(g_flat), 24, 1, "none"))
ap_flat <- e$dance_traj_amp_phase(e$dance_traj_cell_coefs(f_flat, 1))
chk(isTRUE(ap_flat$ok) && any(!ap_flat$table$phase_defined) &&
      all(is.na(ap_flat$table$acrophase_lo[!ap_flat$table$phase_defined])),
    sprintf("with no rhythm, %d of %d cells report NO acrophase interval (Bingham)",
            sum(!ap_flat$table$phase_defined), nrow(ap_flat$table)),
    "an acrophase interval was printed for a cell whose amplitude covers zero")
chk(grepl("UNDEFINED", ap_flat$note %||% ""),
    "and the reason is carried with the result, not left to the reader",
    "the undefined-phase note is missing")

cat("\n-- layer D: phase contrasts wrap ------------------------------------\n")
pc <- e$dance_traj_phase_contrast(co, 1, 2)
chk(!is.null(pc) && abs(pc$diff_time) <= co$effective_period / 2 + 1e-9,
    sprintf("a phase contrast is the WRAPPED difference (%.3f h, within +/- %.1f h)",
            pc$diff_time, co$effective_period / 2),
    "a phase contrast ran outside the half-period it must wrap into")
chk(pc$lo < pc$diff_time && pc$diff_time < pc$hi,
    "and its interval contains it",
    "the contrast interval does not contain the contrast")

cat("\n-- layer B: tau is profiled, and a flat profile says so -------------\n")
# no saturating trend in the data, so the profile MUST be flat: finding A10 made
# quantitative rather than left as a warning beside a number
pr <- e$dance_traj_profile_tau(build(make(41, n_per_group = 6)), 24, 1,
                               c("Group", "Condition"))
chk(isTRUE(pr$ok) && isTRUE(pr$flat) && grepl("NOT identified", pr$message),
    sprintf("with no saturating trend the profile is flat (logLik range %.2f) and says tau is unidentified",
            pr$logLik_range),
    sprintf("a flat tau profile was not detected (range %.2f, flat = %s)",
            pr$logLik_range %||% NA, pr$flat %||% NA))

cat("\n")
if (bad) { cat(sprintf("Trajectory framework tests FAILED  (%d passed, %d failed)\n", ok_n, bad)); quit(status = 1) }
cat(sprintf("Trajectory framework tests PASSED  (%d passed, 0 failed)\n", ok_n))
