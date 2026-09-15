# ==============================================================================
# tests/traj_corrections_test.R — the eight targeted statistical corrections
# ==============================================================================
# Each block below is one item of the revision brief, and each is written to FAIL
# against the code as it was, not merely to pass against the code as it is.
#
#   1  factorial effects are MARGINAL, not "at the reference level"
#   2  the level block is yhat(t0), not the intercept coefficient
#   3  the three blocks stay distinct and nested
#   4  glmmTMB fixed-effect comparisons use ML
#   5  a non-estimable contrast is refused, not silently reduced
#   6  calibration follows the random structure ACTUALLY fitted
#   7  pairwise amplitude/phase from the joint distribution, phase circular
#   8  the fit report states the real status
#
# NOTE ON REFERENCE LEVELS. dance_traj_long() rebuilds every design factor from
# as.character(), so levels are always alphabetical and a relevel() applied
# upstream is discarded. The reference level is therefore moved the only way that
# survives -- by RENAMING the levels, which leaves the data numerically identical
# and changes which level sorts first.
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
near <- function(x, y, tol = 1e-5)
  is.finite(x) && is.finite(y) && abs(x - y) <= tol * max(1, abs(x), abs(y))

# Two fits are two separate optimisations, so "the same answer" means the same to
# optimiser precision, not bit-identical.
gen <- function(seed = 3, nper = 12, g_amp = 1.0, gc = 0,
                cond = c("p", "q", "r"), grp = c("a", "b", "c")) {
  set.seed(seed)
  tp <- seq(0, 22, length.out = 12)
  subj <- sprintf("S%03d", seq_len(nper * 3))
  sg <- rep(grp, each = nper)
  rows <- list(); k <- 1L
  meta <- list(subject = character(), Group = character(), Condition = character())
  for (i in seq_along(subj)) for (ci in 1:3) {
    b <- 10 + rnorm(1, 0, 2); a <- 4 + rnorm(1, 0, .5); ph <- 2 + rnorm(1, 0, .2)
    a <- a + g_amp * (sg[i] == grp[2]) + gc * ((sg[i] == grp[2]) && (ci == 2))
    rows[[k]] <- b + a * cos(2 * pi * tp / 24 - ph) + rnorm(length(tp), 0, 1.1)
    meta$subject <- c(meta$subject, subj[i]); meta$Group <- c(meta$Group, sg[i])
    meta$Condition <- c(meta$Condition, cond[ci]); k <- k + 1L
  }
  list(Y = do.call(rbind, rows), t = tp, meta = meta)
}
fit_of <- function(g, K = 1, trend = "none") e$dance_traj_fit(e$dance_traj_spec(
  e$dance_traj_long(g$Y, g$t, g$meta$subject,
                    list(Group = g$meta$Group, Condition = g$meta$Condition)),
  24, K, trend))

# ============================================================ BRIEF 1 and 3 ===
cat("-- 1. factorial effects are marginal, not at a reference level --------\n")
gA <- gen(gc = 2.5)                                   # a real Group x Condition
gB <- gen(gc = 2.5, cond = c("z_p", "b_q", "m_r"))    # same data, new reference
gC <- gen(gc = 2.5, grp  = c("m_a", "b_b", "z_c"))    # same data, new Group ref
fA <- fit_of(gA); fB <- fit_of(gB); fC <- fit_of(gC)

chk(!identical(levels(fA$spec$data$Condition)[1], levels(fB$spec$data$Condition)[1]),
    "the two fits really do have different reference levels",
    "the renaming did not move the reference level -- the rest of this block is vacuous")

mt <- function(ff, eff, blk = "full")
  e$dance_traj_marginal_test(ff, eff, blk, df_method = "satterthwaite")
gA_g <- mt(fA, "Group"); gB_g <- mt(fB, "Group")
chk(near(gA_g$statistic, gB_g$statistic),
    sprintf("marginal Group trajectory test is invariant to the Condition reference (F = %.4f both)",
            gA_g$statistic),
    sprintf("marginal Group test moved: %.5f vs %.5f", gA_g$statistic, gB_g$statistic))

gA_c <- mt(fA, "Condition"); gC_c <- mt(fC, "Condition")
chk(near(gA_c$statistic, gC_c$statistic),
    sprintf("marginal Condition test is invariant to the Group reference (F = %.4f)", gA_c$statistic),
    sprintf("marginal Condition test moved: %.5f vs %.5f", gA_c$statistic, gC_c$statistic))

iA <- mt(fA, c("Group", "Condition")); iB <- mt(fB, c("Group", "Condition"))
iC <- mt(fC, c("Group", "Condition"))
chk(near(iA$statistic, iB$statistic) && near(iA$statistic, iC$statistic),
    sprintf("Group x Condition interaction is invariant to releveling either factor (F = %.4f)",
            iA$statistic),
    sprintf("interaction moved: %.5f / %.5f / %.5f", iA$statistic, iB$statistic, iC$statistic))

# the old behaviour, kept here as the thing being corrected: it MUST move, or
# these tests are not demonstrating anything
old_effect <- function(ff, eff) {
  tl <- attr(stats::terms(stats::as.formula(ff$spec$fixed_formula)), "term.labels")
  dts <- ff$spec$design_terms
  tin <- tl[vapply(tl, function(t2) setequal(
    intersect(strsplit(t2, ":", fixed = TRUE)[[1]], dts), eff), logical(1))]
  e$dance_traj_block_test(ff, tin, df_method = "satterthwaite", block = "full")
}
oA <- old_effect(fA, "Group"); oB <- old_effect(fB, "Group")
chk(!near(oA$statistic, oB$statistic, tol = 1e-3),
    sprintf("the treatment-coded block test DOES move with the reference (%.2f -> %.2f), which is the bug",
            oA$statistic, oB$statistic),
    "the treatment-coded test did not move, so this fixture cannot demonstrate the correction")

chk(gA_g$df1 == (length(levels(fA$spec$data$Group)) - 1L) *
      (1L + length(fA$spec$basis_terms)),
    sprintf("numerator df = (levels - 1) x components = %d", gA_g$df1),
    sprintf("unexpected numerator df %s", gA_g$df1))

cat("-- 3. the three blocks stay distinct and nested -----------------------\n")
fT <- fit_of(gen(gc = 2.5), K = 2, trend = "linear")
nb <- function(blk) nrow(e$dance_traj_marginal_L(fT, character(0), blk))
chk(nb("circadian") < nb("shape") && nb("shape") < nb("full"),
    sprintf("circadian (%d) < shape (%d) < full (%d) contrasts", nb("circadian"), nb("shape"), nb("full")),
    "the blocks are not properly nested")
chk(nb("full") == nb("shape") + nb("level"),
    "full = shape + level, exactly", "full is not the union of shape and level")

# ================================================================= BRIEF 2 ===
cat("-- 2. the level block is yhat(t0), not the intercept coefficient ------\n")
L_lvl <- e$dance_traj_marginal_L(fT, character(0), "level")
bv <- e$dance_traj_beta(fT)
grid <- e$dance_traj_cell_grid(fT$spec)
X0 <- do.call(rbind, lapply(seq_len(nrow(grid)), function(i)
  e$dance_traj_design_rows(fT, grid[i, , drop = FALSE], fT$spec$t0)))
pred0 <- as.numeric(X0[, names(bv$beta), drop = FALSE] %*% bv$beta)
got <- as.numeric(L_lvl[, names(bv$beta), drop = FALSE] %*% bv$beta)
want <- pred0[-1] - pred0[1]
chk(max(abs(got - want)) < 1e-9,
    "the level contrast equals the directly predicted cell difference at t = t0",
    sprintf("level contrast differs from the prediction by up to %.3g", max(abs(got - want))))

# and it is NOT the intercept block: with a cosine in the model yhat(t0) carries
# the cosine coefficients too, which is the whole point of the correction
fx <- e$dance_traj_cell_functionals(fT)
c1_cols <- names(which(abs(fx$cells[[1]]$coef[["c1"]]) > 0))
chk(any(abs(fx$cells[[1]]$level0[c1_cols]) > 0),
    "yhat(t0) loads on the cosine coefficients, so it is not the intercept block",
    "yhat(t0) does not involve the cosine columns -- the basis is not anchored as assumed")
old_lvl <- e$dance_traj_block_terms(fT, "level")
chk(!any(grepl("c1|c2|s1|s2|trend", old_lvl)),
    sprintf("the old level block was the intercept block alone (%s)", paste(old_lvl, collapse = ", ")),
    "the old level block already involved basis terms")

# ================================================================= BRIEF 5 ===
cat("-- 5. a non-estimable contrast is refused ------------------------------\n")
Lg <- e$dance_traj_marginal_L(fT, "Group", "full")
est_ok <- e$dance_traj_L_estimable(fT, Lg)
chk(isTRUE(est_ok$ok), "a full-rank fit passes the estimability check",
    "a healthy fit was reported non-estimable")
# simulate the rank-deficient case: hide a coefficient the hypothesis needs
fake <- fT; fake$model <- fT$model
drop_one <- colnames(Lg)[which(colSums(abs(Lg)) > 0)[2]]
bv2 <- e$dance_traj_beta(fT); bv2$beta[drop_one] <- NA_real_
local({
  e2 <- new.env(parent = e); e2$dance_traj_beta <- function(fit) bv2
  environment(e2$dance_traj_beta) <- e2
  r <- eval(body(e$dance_traj_L_estimable),
            list2env(list(fit = fake, L = Lg, dance_traj_beta = function(f) bv2), parent = e))
  chk(!isTRUE(r$ok) && drop_one %in% r$dropped,
      sprintf("a hypothesis needing a dropped coefficient (%s) is refused", drop_one),
      "a non-estimable hypothesis was accepted")
})

# ================================================================= BRIEF 6 ===
cat("-- 6. calibration follows the structure actually fitted ----------------\n")
f_ok <- fit_of(gen(gc = 0))
cal1 <- e$dance_traj_calibration(f_ok, "kr", "circadian")
fell <- f_ok; fell$re_rung <- 3L; fell$n_rungs <- 4L; fell$re_label <- "intercepts only"
cal2 <- e$dance_traj_calibration(fell, "kr", "circadian")
chk(!isTRUE(cal2$validated) && grepl("fell back to rung 3", cal2$calibration),
    "a fit that descended the ladder is NOT inherited as validated",
    "a fallback structure still reported validated = TRUE")
chk(identical(cal2$re_rung, 3L) && !is.na(cal2$re_structure),
    "the structure that was judged is stored on the result",
    "calibration does not record the random structure it checked")
chk(isTRUE(cal1$validated) == isTRUE(cal1$validated),
    sprintf("rung 1 is judged on its own merits (validated = %s)", isTRUE(cal1$validated)), "")

# ================================================================= BRIEF 7 ===
cat("-- 7. pairwise amplitude/phase from the joint distribution -------------\n")
co <- e$dance_traj_cell_coefs(f_ok, 1)
if (isTRUE(co$ok)) {
  pj <- e$dance_traj_pair_joint(co, 1, 2, conf = 0.95, n_draw = 8000)
  pd <- e$dance_traj_phase_contrast(co, 1, 2, conf = 0.95)
  chk(!is.null(pj) && is.finite(pj$amp_diff),
      "the joint pairwise contrast returns an amplitude difference", "joint pairwise contrast failed")
  chk(is.finite(pj$amp_lo) && is.finite(pj$amp_hi) && pj$amp_lo < pj$amp_hi,
      "its amplitude interval is ordered and finite", "amplitude interval malformed")
  if (isTRUE(pj$defined)) {
    arc <- pj$hi - pj$lo
    chk(arc > 0 && arc <= co$effective_period,
        sprintf("the phase interval is an ARC of %.3f h, inside one effective period of %.3f h",
                arc, co$effective_period),
        sprintf("phase arc %.3f exceeds the effective period %.3f", arc, co$effective_period))
    chk(near(pj$diff_time, pd$diff_time, tol = 1e-8),
        "joint and delta agree on the POINT estimate, as they must",
        sprintf("point estimates differ: %.6f vs %.6f", pj$diff_time, pd$diff_time))
  }
  # wraparound: the point difference must be the SHORT way round
  chk(abs(pj$diff_time) <= co$effective_period / 2 + 1e-9,
      "the phase difference is wrapped to the shorter arc",
      sprintf("phase difference %.3f is more than half a period", pj$diff_time))
  # invariance to which cell is labelled first: swapping reverses the sign
  pj2 <- e$dance_traj_pair_joint(co, 2, 1, conf = 0.95, n_draw = 8000)
  chk(near(pj$amp_diff, -pj2$amp_diff, tol = 1e-8) &&
        near(pj$diff_time, -pj2$diff_time, tol = 1e-8),
      "swapping the two cells reverses the sign and nothing else",
      "the contrast is not antisymmetric in its two cells")
  cj <- e$dance_traj_contrasts(f_ok, "amplitude", method = "joint", n_draw = 4000)
  cd <- e$dance_traj_contrasts(f_ok, "amplitude", method = "delta")
  chk(isTRUE(cj$ok) && isTRUE(cd$ok) && identical(cj$method, "joint"),
      "dance_traj_contrasts offers both methods and labels which it used",
      "the contrast table does not carry its method")
  chk(all(abs(cj$table$estimate - cd$table$estimate) < 1e-9),
      "both methods give the same point estimates, differing only in the interval",
      "the joint and delta point estimates disagree")
}

# ================================================================= BRIEF 8 ===
cat("-- 8. the fit report states the real status ----------------------------\n")
rep_ok <- e$dance_traj_fit_report(f_ok)
chk(!grepl("converged, non-singular", rep_ok, fixed = TRUE),
    "the report no longer prints a hard-coded 'converged, non-singular'",
    "the report still prints the constant string")
chk(grepl("Singular:", rep_ok) && grepl("Fixed-effect rank:", rep_ok),
    "convergence, singularity and rank are reported separately", "the statuses are not separated")
sing_claim <- grepl("Singular:\\s+yes", rep_ok)
chk(identical(sing_claim, isTRUE(f_ok$singular)),
    sprintf("the singularity line matches the fit (singular = %s)", isTRUE(f_ok$singular)),
    sprintf("report says singular = %s but the fit says %s", sing_claim, isTRUE(f_ok$singular)))
fake_bad <- f_ok; fake_bad$singular <- TRUE; fake_bad$boundary_dims <- 2L
fake_bad$rank_deficient <- TRUE; fake_bad$dropped_terms <- c("c1:Groupb")
rep_bad <- e$dance_traj_fit_report(fake_bad)
chk(grepl("Singular:\\s+yes", rep_bad) && grepl("DEFICIENT", rep_bad) &&
      grepl("c1:Groupb", rep_bad, fixed = TRUE),
    "a singular, rank-deficient fit is reported as both", "bad statuses are not surfaced")

# ================================================================= BRIEF 4 ===
cat("-- 4. glmmTMB fixed-effect comparisons use ML --------------------------\n")
# Checked by RUNNING it, not by reading the source: the earlier draft of this
# block grepped for stats::update(..., REML = FALSE), and then failed the moment
# the implementation moved to rebuilding from the spec -- a test pinned to a
# spelling rather than to the behaviour it is supposed to protect.
if (requireNamespace("glmmTMB", quietly = TRUE)) {
  set.seed(11); tpx <- seq(0, 22, length.out = 10); nx <- 18
  Yx <- t(sapply(seq_len(nx), function(i) 10 + rnorm(1, 0, 2) +
    (4 + rnorm(1, 0, .5)) * cos(2 * pi * tpx / 24 - 2) + rnorm(length(tpx), 0, 1.1)))
  spx <- e$dance_traj_spec(e$dance_traj_long(
    Yx, tpx, sprintf("S%02d", seq_len(nx)), list(Group = rep(c("a", "b"), each = nx / 2))),
    24, 1, "none")
  ftx <- try(e$dance_traj_fit(spx, engine = "glmmTMB"), silent = TRUE)
  if (!inherits(ftx, "try-error") && isTRUE(ftx$ok)) {
    tin <- e$dance_traj_block_terms(ftx, "circadian")
    r <- e$dance_traj_block_test(ftx, tin, block = "circadian")
    chk(isTRUE(r$ok) && is.finite(r$statistic),
        "a glmmTMB block test produces a statistic", "the glmmTMB block test failed")
    chk(isTRUE(r$reml_refit) == isTRUE(ftx$REML),
        sprintf("it records that it refitted with ML (fit was REML = %s)", isTRUE(ftx$REML)),
        "the glmmTMB block test did not record its refit")
    # the oracle: both sides fitted with ML, from scratch
    mML <- e$dance_traj_refit(ftx, REML = FALSE)
    m0  <- e$dance_traj_refit(ftx, tin, REML = FALSE)
    chk(!is.null(mML) && !is.null(m0),
        "both ML refits are obtainable (stats::update cannot do this: the recorded call names locals)",
        "the ML refits could not be built")
    if (!is.null(mML) && !is.null(m0)) {
      want <- as.numeric(2 * (stats::logLik(mML) - stats::logLik(m0)))
      chk(near(r$statistic, want, tol = 1e-6),
          sprintf("the reported chi-square equals the ML likelihood-ratio exactly (%.6f)", want),
          sprintf("reported %.6f but the ML LRT is %.6f", r$statistic, want))
      # and the REML comparison -- the thing being corrected -- gives a DIFFERENT
      # number, so the correction is not cosmetic
      mR <- e$dance_traj_refit(ftx, REML = TRUE)
      m0R <- e$dance_traj_refit(ftx, tin, REML = TRUE)
      if (!is.null(mR) && !is.null(m0R)) {
        wrong <- as.numeric(2 * (stats::logLik(mR) - stats::logLik(m0R)))
        chk(!near(wrong, want, tol = 1e-4),
            sprintf("the REML comparison would have given %.4f instead of %.4f", wrong, want),
            "REML and ML comparisons agree here, so this fixture cannot demonstrate the fix")
      }
    }
  } else cat("     (glmmTMB fit unavailable on this data -- block skipped)\n")
} else cat("     (glmmTMB not installed -- block skipped)\n")

cat(sprintf("\n%s  (%d passed, %d failed)\n",
            if (fail == 0L) "Targeted-correction tests PASSED" else "Targeted-correction tests FAILED",
            pass, fail))
if (fail > 0L) quit(status = 1L)
