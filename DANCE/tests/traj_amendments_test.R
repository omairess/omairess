# ==============================================================================
# tests/traj_amendments_test.R — the P21 phase-3 statistical amendments
# ==============================================================================
# One check per amendment that can be checked, plus the refusals. What is NOT
# here is the expanded calibration grid: that is the pre-phase-4 gate and needs
# thousands of simulations, not a unit test.
# ==============================================================================

`%||%` <- function(a, b) if (is.null(a)) b else a
app <- if (dir.exists("server")) "." else ".."
e <- new.env()
for (f in c("08_helpers_cosinor", "08c_helpers_circstat", "08d_helpers_traj",
            "08e_helpers_trajfit", "08f_helpers_trajinf", "08g_helpers_trajband"))
  sys.source(file.path(app, "server", paste0(f, ".R")), envir = e)

pass <- 0L; fail <- 0L
chk <- function(ok, a, b) { if (isTRUE(ok)) { pass <<- pass + 1L; cat("ok   ", a, "\n") }
                            else { fail <<- fail + 1L; cat("FAIL ", b, "\n") } }

# a generator with an arbitrary number of within-factor levels
gen <- function(seed, n_per_group = 6, groups = c("a", "b"), conds = c("p", "q"),
                visits = NULL, nt = 12, amp = 4, noise = 1, times = NULL) {
  set.seed(seed)
  tp <- times %||% seq(0, 22, length.out = nt)
  subj <- sprintf("S%03d", seq_len(n_per_group * length(groups)))
  sg <- rep(groups, each = n_per_group)
  rows <- list(); k <- 1L
  meta <- list(subject = character(0), Group = character(0),
               Condition = character(0), Visit = character(0))
  vs <- visits %||% "v1"
  for (i in seq_along(subj)) {
    sb <- rnorm(1, 0, 2); sa <- rnorm(1, 0, .4); sp <- rnorm(1, 0, .2)
    for (cc in conds) for (vv in vs) {
      base <- 10 + sb + rnorm(1, 0, 2); a <- amp + sa + rnorm(1, 0, .4)
      ph <- 2 + sp + rnorm(1, 0, .2)
      rows[[k]] <- base + a * cos(2 * pi * tp / 24 - ph) + rnorm(length(tp), 0, noise)
      meta$subject <- c(meta$subject, subj[i]); meta$Group <- c(meta$Group, sg[i])
      meta$Condition <- c(meta$Condition, cc); meta$Visit <- c(meta$Visit, vv)
      k <- k + 1L
    }
  }
  list(Y = do.call(rbind, rows), t = tp, meta = meta)
}
build <- function(g, with_visit = FALSE) {
  fl <- list(Group = g$meta$Group, Condition = g$meta$Condition)
  if (with_visit) fl$Visit <- g$meta$Visit
  e$dance_traj_long(g$Y, g$t, g$meta$subject, fl)
}

suppressWarnings(suppressMessages({

cat("-- A2: the curve id is general, not subject:Condition ----------------\n")
d2 <- build(gen(1))
sp2 <- e$dance_traj_spec(d2, 24, 1, "none")
chk(isTRUE(sp2$ok) && sp2$n_curves == 24 && sp2$n_participants == 12 &&
      grepl("participant x Condition", sp2$curve_definition),
    sprintf("2 conditions x 12 participants = %d curves, defined as %s",
            sp2$n_curves, sp2$curve_definition),
    "the curve count or definition is wrong for a 2-condition design")
d3 <- build(gen(2, conds = c("p", "q", "r")))
sp3 <- e$dance_traj_spec(d3, 24, 1, "none")
chk(isTRUE(sp3$ok) && sp3$n_curves == 36,
    sprintf("a THREE-level within factor gives %d curves from the same code", sp3$n_curves),
    "a 3-level within factor did not generalise")
d4 <- build(gen(3, conds = c("p", "q"), visits = c("v1", "v2")), with_visit = TRUE)
sp4 <- e$dance_traj_spec(d4, 24, 1, "none")
chk(isTRUE(sp4$ok) && sp4$n_curves == 48 &&
      grepl("Condition", sp4$curve_definition) && grepl("Visit", sp4$curve_definition),
    sprintf("TWO within factors compose into one curve column (%d curves: %s)",
            sp4$n_curves, sp4$curve_definition),
    "two within factors did not compose into one curve level")
chk(grepl("| curve)", sp4$re_ladder[[1]]$formula, fixed = TRUE) &&
      !grepl("subject:", sp4$re_ladder[[1]]$formula, fixed = TRUE),
    "and the ladder groups by that column, never by a spelled-out interaction",
    sprintf("the ladder still spells out an interaction: %s", sp4$re_ladder[[1]]$formula))
# a purely between design has no curve level to speak of
db <- e$dance_traj_long(gen(4, conds = "p")$Y, gen(4, conds = "p")$t,
                        gen(4, conds = "p")$meta$subject,
                        list(Group = gen(4, conds = "p")$meta$Group))
spb <- e$dance_traj_spec(db, 24, 1, "none")
chk(isTRUE(spb$ok) && isFALSE(spb$has_curve_level) &&
      spb$n_curves == spb$n_participants,
    "with no within factor a curve IS a participant, and no curve rung is offered",
    "a between-only design was given a curve level")

cat("\n-- a single-level factor is not a factor ------------------------------\n")
# Found by the validation grid, not by review: a within-only design (one group)
# produced NO FITS AT ALL. A one-level factor satisfies "every participant sees
# exactly one level" and so read as between-participant, and then every column
# it contributed was aliased with the intercept. This happens in practice
# whenever a filter leaves one group standing.
g1 <- gen(8, n_per_group = 16, groups = "a", conds = c("p", "q", "r"))
d1 <- build(g1)
c1 <- e$dance_traj_classify(d1)
chk(c1$role[c1$factor == "Group"] == "constant" && c1$n_levels[c1$factor == "Group"] == 1,
    "a one-level factor is classified 'constant', not 'between'",
    sprintf("a one-level factor was classified '%s'", c1$role[c1$factor == "Group"]))
sp1 <- e$dance_traj_spec(d1, 24, 1, "none")
chk(isTRUE(sp1$ok) && !grepl("Group", sp1$fixed_formula) &&
      identical(sp1$constant_terms, "Group"),
    sprintf("it is dropped from the formula (%s) and the drop is RECORDED", sp1$fixed_formula),
    "a constant factor survived into the model, or its removal was silent")
chk(grepl("aliased with the intercept", sp1$constant_note %||% "") &&
      grepl("Dropped:", e$dance_traj_describe(sp1), fixed = TRUE),
    "and the readout says which factor went and why",
    "the constant-factor drop is not explained in the readout")
f1 <- e$dance_traj_fit(sp1)
o1 <- e$dance_traj_omnibus(f1, "full")
chk(isTRUE(f1$ok) && isTRUE(o1$ok) && is.finite(o1$p),
    sprintf("the within-only design now fits and its full block is estimable (p = %.3f)", o1$p),
    "a within-participant-only design still produces no usable test")

cat("\n-- A10: the classification is shown, and can be overridden -----------\n")
desc <- e$dance_traj_describe(sp2)
chk(grepl("Group", desc) && grepl("between", desc) && grepl("read from the data", desc) &&
      grepl("Curve (session)", desc, fixed = TRUE),
    "the readout lists every factor's role, where it came from, and the curve definition",
    "the description does not expose the classification")
sp_o <- e$dance_traj_spec(d2, 24, 1, "none", roles = c(Condition = "between"))
chk(isTRUE(sp_o$ok) && sp_o$classification$role[sp_o$classification$factor == "Condition"] == "between" &&
      sp_o$classification$source[sp_o$classification$factor == "Condition"] == "user" &&
      isFALSE(sp_o$has_curve_level),
    "a user override changes the role, is recorded as an override, and changes the ladder",
    "the role override did not take effect or was not recorded")
chk(grepl("SET BY YOU", e$dance_traj_describe(sp_o), fixed = TRUE),
    "and the readout says the role was set by the user, not read from the data",
    "an overridden role is not flagged in the readout")
chk(!isTRUE(e$dance_traj_spec(d2, 24, 1, "none", roles = c(Condition = "nonsense"))$ok),
    "an unknown role is refused rather than silently ignored",
    "an invalid role was accepted")

cat("\n-- A1: three omnibus questions, nested -------------------------------\n")
fit <- e$dance_traj_fit(sp2)
tf <- e$dance_traj_block_terms(fit, "full")
ts <- e$dance_traj_block_terms(fit, "shape")
tc <- e$dance_traj_block_terms(fit, "circadian")
tl <- e$dance_traj_block_terms(fit, "level")
chk(all(tc %in% ts) && all(ts %in% tf) && all(tl %in% tf) &&
      length(intersect(tl, ts)) == 0,
    sprintf("full (%d) > shape (%d) > circadian (%d); level (%d) is inside full and disjoint from shape",
            length(tf), length(ts), length(tc), length(tl)),
    "the block definitions are not nested as amendment 1 specifies")
chk(setequal(tf, union(ts, tl)),
    "and full is exactly shape plus level -- nothing falls between the definitions",
    "full is not the union of shape and level, so some term belongs to no question")

cat("\n-- A5: no N = 50 statistical cutoff ----------------------------------\n")
o <- e$dance_traj_omnibus(fit, "circadian")
chk(grepl("Kenward-Roger", o$method) && !grepl("<=", o$method) && !grepl("50", o$method),
    sprintf("the method string names the test without a threshold: '%s'", o$method),
    sprintf("the method string still quotes a cutoff: '%s'", o$method))
os <- e$dance_traj_omnibus(fit, "circadian", df_method = "satterthwaite")
chk(grepl("Satterthwaite", os$method) && grepl("cheaper", os$method_note %||% ""),
    "Satterthwaite can be requested explicitly and says it is the cheaper alternative",
    "requesting Satterthwaite did not work or was not labelled")
oc <- e$dance_traj_omnibus(fit, "circadian", kr_max_subjects = 4)
chk(grepl("Satterthwaite", oc$method) && grepl("COMPUTATIONAL", oc$method_note %||% ""),
    "a runtime cap is honoured but labelled a COMPUTATIONAL policy, not a statistical rule",
    "the runtime cap is not labelled as a computational policy")

cat("\n-- A4: the calibration does not propagate ----------------------------\n")
chk(isTRUE(o$validated) && isTRUE(o$provisional) &&
      grepl("PROVISIONALLY", o$calibration) && grepl("[0.011, 0.099]", o$calibration, fixed = TRUE),
    "the circadian block calls itself PROVISIONAL and quotes exact binomial intervals",
    "the calibration statement is missing its provisional status or its intervals")
chk(isFALSE(e$dance_traj_omnibus(fit, "full")$validated) &&
      grepl("only the circadian block", e$dance_traj_omnibus(fit, "full")$calibration),
    "the FULL block does not borrow the circadian block's evidence",
    "an unsimulated block claimed the circadian block's calibration")
chk(isFALSE(os$validated) && grepl("Satterthwaite", os$calibration),
    "and neither does a Satterthwaite fit",
    "a Satterthwaite result claimed the Kenward-Roger grid's calibration")

cat("\n-- A3: five statuses, kept apart -------------------------------------\n")
st <- e$dance_traj_status(fit$model)
chk(all(c("converged", "singular", "boundary_dims", "rank_deficient",
          "optimizer_failure") %in% names(st)),
    "converged, singular, boundary, rank-deficient and optimizer failure are separate fields",
    "the statuses are still collapsed into one")
chk(is.list(fit$re_collapse) && all(vapply(fit$re_collapse, function(x)
      all(c("n_dim", "n_collapsed", "sdev") %in% names(x)), logical(1))),
    sprintf("rePCA names which dimensions collapsed (%d of %d at the participant level)",
            fit$re_collapse[[1]]$n_collapsed, fit$re_collapse[[1]]$n_dim),
    "the collapsed random-effect dimensions are not reported")
chk(grepl("NOT a claim", fit$note %||% "not singular here", fixed = TRUE) ||
      !isTRUE(fit$singular),
    "a retained singular fit says so WITHOUT claiming the maximal structure is right",
    "a singular fit was reported as though retention settled the question")

cat("\n-- A7: discrete AR(1) is refused on an irregular grid -----------------\n")
d_irr <- build(gen(5, times = c(0, 1, 2, 4, 8, 12, 13, 14, 18, 22)))
sp_irr <- e$dance_traj_spec(d_irr, 24, 1, "none")
chk(isTRUE(sp_irr$ok) && isFALSE(sp_irr$time_regular),
    sprintf("the spec detects uneven sampling (gaps %.3g to %.3g)",
            sp_irr$time_gaps[1], sp_irr$time_gaps[2]),
    "uneven sampling was not detected")
r_ar1 <- e$dance_traj_fit(sp_irr, engine = "glmmTMB", residual_cor = "ar1")
chk(isFALSE(r_ar1$ok) && grepl("NOT evenly", r_ar1$message) && grepl("'ou'", r_ar1$message),
    "a discrete AR(1) on that grid is REFUSED, naming the continuous-time alternative",
    "a discrete AR(1) was accepted on an irregular grid")
chk(isTRUE(e$dance_traj_spec(d2, 24, 1, "none")$time_regular),
    "and an even grid is still recognised as even",
    "an even grid was misread as uneven")
r_lmer <- e$dance_traj_fit(sp2, engine = "lmer", residual_cor = "ou")
chk(isFALSE(r_lmer$ok) && grepl("glmmTMB", r_lmer$message),
    "lme4 refuses any residual correlation structure and names the engine that can",
    "lme4 accepted a residual correlation structure")

cat("\n-- A6: joint-distribution amplitude/phase ----------------------------\n")
co <- e$dance_traj_cell_coefs(fit, 1)
aj <- e$dance_traj_amp_phase_joint(co, n_draw = 8000)
ad <- e$dance_traj_amp_phase(co)
chk(isTRUE(aj$ok) && all(aj$table$amplitude_lo > 0),
    "the joint method's amplitude interval cannot go below zero",
    "a joint amplitude interval reached zero or below")
chk(max(abs(aj$table$amplitude - ad$table$amplitude)) < 1e-12,
    "both methods agree exactly on the point estimate",
    "the two methods disagree on the point estimate")
chk(max(abs(aj$table$amplitude_lo - ad$table$amplitude_lo)) < 0.3,
    sprintf("and on the interval, where the vector is far from zero (max diff %.3f)",
            max(abs(aj$table$amplitude_lo - ad$table$amplitude_lo))),
    "the two methods disagree badly in the regime where the delta method should be fine")
aj2 <- e$dance_traj_amp_phase_joint(co, n_draw = 8000)
chk(identical(aj$table$amplitude_lo, aj2$table$amplitude_lo),
    "the simulated interval is reproducible: the same fit gives the same number twice",
    "the simulated interval moved between two runs of the same fit")
# near the origin, where the delta method is wrong and the joint method must say so
fit0 <- e$dance_traj_fit(e$dance_traj_spec(build(gen(6, amp = 0, noise = 3)), 24, 1, "none"))
co0 <- e$dance_traj_cell_coefs(fit0, 1)
aj0 <- e$dance_traj_amp_phase_joint(co0, n_draw = 8000)
chk(isTRUE(aj0$ok) && any(!aj0$table$phase_defined) &&
      all(is.na(aj0$table$acrophase_lo[!aj0$table$phase_defined])),
    sprintf("with no rhythm, %d of %d cells report NO phase interval",
            sum(!aj0$table$phase_defined), nrow(aj0$table)),
    "a phase interval was reported for a vector that may be the zero vector")
chk(grepl("ARC", aj$note, fixed = TRUE) && "acrophase_arc_time" %in% names(aj$table),
    "the phase interval is reported as an arc with its length, not as a symmetric interval",
    "the phase interval is not labelled as an arc")

cat("\n-- A11: simple effects and interaction contrasts ---------------------\n")
se <- e$dance_traj_simple_effects(fit, "Condition", at = list(Group = "a"))
chk(isTRUE(se$ok) && nrow(se$table) == 1 &&
      all(c(se$table$cell1, se$table$cell2) %in% c("a x p", "a x q")),
    "the simple effect of Condition at Group = a is the single contrast inside that slice",
    "the simple effect did not restrict to its slice")
chk(grepl("not by refitting", se$note) && grepl("this slice only", se$note),
    "and it says it came from the one fitted model, adjusted within the slice only",
    "the simple effect does not state its provenance or its family")
chk(isFALSE(e$dance_traj_simple_effects(fit, "Nonesuch", at = list(Group = "a"))$ok) &&
      isFALSE(e$dance_traj_simple_effects(fit, "Condition", at = list())$ok),
    "an unknown factor, and an incomplete `at`, are both refused",
    "a malformed simple-effect request was accepted")
ic <- e$dance_traj_interaction_contrast(fit, "a x p", "a x q", "b x p", "b x q")
chk(isTRUE(ic$ok) && is.finite(ic$se) && grepl("difference of differences", ic$note),
    sprintf("the interaction contrast is one linear combination (estimate %.3f, se %.3f)",
            ic$estimate, ic$se),
    "the interaction contrast failed")
chk(grepl("NOT answered by observing that one difference is significant", ic$note, fixed = TRUE),
    "and it says plainly what an interaction is not",
    "the interaction note does not warn against the two-p-values fallacy")

cat("\n-- A12: bootstrap copies get NEW cluster ids -------------------------\n")
subj <- rep(c("s1", "s2", "s3"), each = 4)
cond <- rep(c("p", "p", "q", "q"), 3)
set.seed(9)
bc <- e$dance_boot_clusters(subj, cond, index = c(1L, 1L, 1L))
chk(length(unique(bc$subject)) == 3 && length(bc$rows) == 12,
    "participant 1 drawn three times becomes THREE clusters, not one with tripled rows",
    sprintf("three draws of one participant gave %d cluster(s)", length(unique(bc$subject))))
chk(all(table(bc$subject) == 4) &&
      all(vapply(split(cond[bc$rows], bc$subject), function(x)
            setequal(x, c("p", "q")), logical(1))),
    "and each copy keeps that participant's complete set of conditions",
    "a bootstrap copy lost part of its participant's within-subject set")
chk(length(unique(bc$curve)) == 6,
    "the curve ids are rebuilt from the new participant ids, so copies do not collide",
    "the bootstrap curve ids collide across copies")

cat("\n-- A13: participant curves are labelled as shrunken ------------------\n")
pc <- e$dance_traj_participant_curves(fit, 1)
if (isTRUE(pc$ok)) {
  chk(isTRUE(pc$shrunken) && grepl("shrunken", pc$label) &&
        grepl("UNDERSTATES", pc$note),
      sprintf("BLUP amplitudes are labelled shrunken; their SD is %.3f against a model SD of %.3f",
              pc$sd_blup_amplitude, pc$sd_population_amplitude),
      "participant curves came back without the shrinkage warning")
} else {
  chk(grepl("no participant-specific harmonic terms", pc$message),
      "with no participant-level harmonic terms it refuses rather than inventing curves",
      "participant curves failed for the wrong reason")
}

cat("\n-- A9: selection is separated from confirmation ----------------------\n")
pre <- e$dance_traj_selection_state()
chk(identical(pre$mode, "pre-specified") && isFALSE(pre$post_selection),
    "a pre-specified model says so",
    "the pre-specified state is wrong")
post <- e$dance_traj_selection_state(c("harmonics", "trend"))
chk(identical(post$mode, "data-selected") && isTRUE(post$post_selection) &&
      grepl("CONDITIONAL ON THE SELECTED MODEL", post$note) &&
      grepl("No correction is applied", post$note),
    "a data-selected model says the p-values are conditional, and that nothing was corrected",
    "the data-selected state does not warn about post-selection inference")
sf <- e$dance_traj_selection_from_fit(fit)
chk(identical(sf$post_selection, isTRUE(fit$simplified)),
    sprintf("and descending the ladder counts as selection (this fit: rung %d, %s)",
            fit$re_rung, sf$mode),
    "a ladder simplification was not counted as model selection")

cat("\n-- A8: tau uncertainty does not vanish -------------------------------\n")
pt <- e$dance_traj_profile_tau(build(gen(7)), 24, 1, c("Group", "Condition"),
                               tau_grid = exp(seq(log(2), log(40), length.out = 6)))
chk(isTRUE(pt$ok) && isTRUE(pt$tau_estimated),
    "the profile records that tau was estimated",
    "the tau profile does not record that tau was estimated")
if (isTRUE(pt$flat)) {
  chk(is.null(pt$fit) && grepl("NOT identified", pt$message),
      "a flat profile SUPPRESSES the fit rather than handing back a point on the ridge",
      "a flat profile still returned a fitted model")
} else {
  chk(grepl("CONDITIONAL ON THE SELECTED VALUE", pt$fit$tau_warning %||% ""),
      "and a fit at the selected tau carries the conditional-inference warning",
      "the fit at the selected tau does not warn that tau was estimated")
}

}))

cat(sprintf("\n%s  (%d passed, %d failed)\n",
            if (fail == 0) "Phase-3 amendment tests PASSED" else "FAILURES", pass, fail))
if (fail > 0) quit(status = 1)
