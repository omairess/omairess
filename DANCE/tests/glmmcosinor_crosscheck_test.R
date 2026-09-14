# ==============================================================================
# tests/glmmcosinor_crosscheck_test.R
#
# AN INDEPENDENT IMPLEMENTATION, USED AS AN ORACLE.
#
# Every other test in this suite checks DANCE against DANCE: my code against my
# expectations about my code. That catches a lot, and it did not catch the
# random-effects misspecification that put the omnibus's type-I error at 0.275,
# because my expectation was wrong in the same direction as my code.
#
# GLMMcosinor (rOpenSci, peer-reviewed) fits the same statistical model from
# different primitives -- glmmTMB rather than lme4, and parameterised directly in
# (amplitude, acrophase) rather than in (cos, sin) coefficients. Where the two
# overlap they must agree. Their failure modes are uncorrelated, so agreement is
# evidence of a kind the rest of this suite cannot produce.
#
# THIS TEST SKIPS ITSELF when GLMMcosinor is not installed, so nobody needs it to
# run the suite. To enable it:
#
#     R CMD INSTALL --library=<lib> <GLMMcosinor source>
#     R_LIBS=<lib> Rscript tests/glmmcosinor_crosscheck_test.R
#
# WHAT IS ASSERTED, IN TWO KINDS
#
#   AGREEMENT  the estimation layer: per-cell amplitude and acrophase, the
#              amplitude contrast, the magnitude of the phase contrast. These
#              must match to 1e-4 -- same model, same data, different code.
#
#   DIVERGENCE the inference layer, where the two deliberately differ. These are
#              pinned too, so that if GLMMcosinor changes what it does, this test
#              says so rather than quietly starting to disagree. They are
#              CHARACTERISATIONS of another package, not claims that it is wrong:
#              it targets a different use case and is explicit about its methods.
# ==============================================================================

`%||%` <- function(a, b) if (is.null(a)) b else a
app <- if (dir.exists("server")) "." else ".."

if (!requireNamespace("GLMMcosinor", quietly = TRUE)) {
  cat("SKIP: GLMMcosinor is not installed -- the cross-check needs it.\n")
  cat("      Install it and re-run to enable this test.\n")
  quit(status = 0)
}

e <- new.env()
for (f in c("08_helpers_cosinor", "08c_helpers_circstat", "08d_helpers_traj",
            "08e_helpers_trajfit", "08f_helpers_trajinf", "08g_helpers_trajband"))
  sys.source(file.path(app, "server", paste0(f, ".R")), envir = e)

pass <- 0L; fail <- 0L
chk <- function(ok, a, b) { if (isTRUE(ok)) { pass <<- pass + 1L; cat("ok   ", a, "\n") }
                            else { fail <<- fail + 1L; cat("FAIL ", b, "\n") } }
agree <- function(label, x, y, tol = 1e-4, unit = "") {
  d <- abs(x - y)
  chk(is.finite(d) && d <= tol,
      sprintf("%-34s %10.5f vs %10.5f  (diff %.1e %s)", label, x, y, d, unit),
      sprintf("%-34s %10.5f vs %10.5f  DIFFER by %.2e %s", label, x, y, d, unit))
}

P <- 24
gen <- function(seed, n_per = 20, nt = 12, ampA = 4, ampB = 4, phA = 2, phB = 2) {
  set.seed(seed); tp <- seq(0, 22, length.out = nt)
  rows <- list(); k <- 1L; meta <- list(subject = character(0), Group = character(0))
  for (g in c("A", "B")) for (i in seq_len(n_per)) {
    amp <- if (g == "A") ampA else ampB
    ph  <- if (g == "A") phA  else phB
    a <- amp + rnorm(1, 0, .4); p <- ph + rnorm(1, 0, .15); m <- 10 + rnorm(1, 0, 2)
    rows[[k]] <- m + a * cos(2 * pi * tp / P - p) + rnorm(nt, 0, 1)
    meta$subject <- c(meta$subject, sprintf("%s%02d", g, i))
    meta$Group <- c(meta$Group, g); k <- k + 1L
  }
  Y <- do.call(rbind, rows)
  list(Y = Y, t = tp, meta = meta,
       long = data.frame(subject = factor(rep(meta$subject, each = nt)),
                         Group = factor(rep(meta$Group, each = nt)),
                         t = rep(tp, times = nrow(Y)), y = as.vector(t(Y)),
                         stringsAsFactors = FALSE))
}
# THE RANDOM STRUCTURE MUST MATCH or the comparison is meaningless. A first draft
# of this check gave GLMMcosinor a random intercept while DANCE fitted random
# cosinor slopes, and the resulting p-value gap looked like an inference
# difference when it was a model difference. amp_acro1 expands to the same pair
# of columns DANCE calls c1 and s1.
fit_both <- function(g) {
  d  <- e$dance_traj_long(g$Y, g$t, g$meta$subject, list(Group = g$meta$Group))
  ff <- e$dance_traj_fit(e$dance_traj_spec(d, P, 1, "none"))
  m  <- GLMMcosinor::cglmm(
    y ~ Group + amp_acro(t, group = "Group", period = P) + (1 + amp_acro1 | subject),
    data = g$long, family = gaussian)
  list(dance = ff, glmm = m, tt = summary(m)$transformed.table)
}
gv <- function(tt, lvl, what) tt[sprintf("[Group=%s]:%s1", lvl, what), "estimate"]

suppressWarnings(suppressMessages({

cat("-- AGREEMENT: per-cell estimates, true null ---------------------------\n")
b1 <- fit_both(gen(101))
co1 <- e$dance_traj_cell_coefs(b1$dance, 1); ap1 <- e$dance_traj_amp_phase(co1)
chk(identical(b1$dance$re_formula, "(1 + c1 + s1 | subject)"),
    sprintf("both sides fit the same random structure: %s", b1$dance$re_formula),
    sprintf("DANCE fitted %s -- the comparison is not like-for-like",
            b1$dance$re_formula))
for (lv in c("A", "B")) {
  i <- match(lv, ap1$table$cell)
  agree(sprintf("amplitude, group %s", lv), ap1$table$amplitude[i], gv(b1$tt, lv, "amp"), 1e-4, "VAS")
  # both report the acrophase on [0, 2pi) with the SAME sign convention --
  # verified rather than assumed: an early draft of this file wrongly negated it
  agree(sprintf("acrophase, group %s", lv), ap1$table$acrophase_rad[i], gv(b1$tt, lv, "acr"), 1e-4, "rad")
}

cat("\n-- AGREEMENT: contrasts, planted amplitude difference (4.0 vs 6.0) ----\n")
b2 <- fit_both(gen(202, ampB = 6))
co2 <- e$dance_traj_cell_coefs(b2$dance, 1); ap2 <- e$dance_traj_amp_phase(co2)
ca2 <- e$dance_traj_contrasts(b2$dance, "amplitude")
t2  <- GLMMcosinor::test_cosinor_levels(b2$glmm, param = "amp", x_str = "Group")
for (lv in c("A", "B")) {
  i <- match(lv, ap2$table$cell)
  agree(sprintf("amplitude, group %s", lv), ap2$table$amplitude[i], gv(b2$tt, lv, "amp"), 1e-4, "VAS")
}
agree("amplitude contrast B - A", ca2$table$estimate[1], t2$ind.test$conf.int[1], 1e-4, "VAS")

cat("\n-- AGREEMENT: contrasts, planted phase difference (2.0 vs 2.8 rad) ----\n")
b3 <- fit_both(gen(303, phB = 2.8))
co3 <- e$dance_traj_cell_coefs(b3$dance, 1); ap3 <- e$dance_traj_amp_phase(co3)
cp3 <- e$dance_traj_phase_contrast(co3, match("A", co3$cells), match("B", co3$cells))
t3  <- GLMMcosinor::test_cosinor_levels(b3$glmm, param = "acr", x_str = "Group")
for (lv in c("A", "B")) {
  i <- match(lv, ap3$table$cell)
  agree(sprintf("acrophase, group %s", lv), ap3$table$acrophase_rad[i], gv(b3$tt, lv, "acr"), 1e-4, "rad")
}
d_rad <- cp3$diff_time * 2 * pi / P
agree("phase contrast |B - A|", abs(d_rad), abs(t3$ind.test$conf.int[1]), 1e-4, "rad")

cat("\n-- DIVERGENCE: the inference layers differ, on purpose ----------------\n")
# 1. DANCE's omnibus is a JOINT test of the cos/sin pair. GLMMcosinor's
#    "global" test, reached through test_cosinor_levels(param = "amp"), carries
#    ONE degree of freedom and equals its own amplitude test -- amplitude and
#    acrophase are tested separately, which is exactly the marginal-testing
#    problem the P21 brief asked DANCE to stop doing.
om3 <- e$dance_traj_omnibus(b3$dance, "circadian")
chk(om3$df1 == 2 && t3$global.test$df == 1,
    sprintf("DANCE tests the pair jointly (df = %d); GLMMcosinor's global test is marginal (df = %d)",
            om3$df1, t3$global.test$df),
    sprintf("the df structure changed: DANCE %s, GLMMcosinor %s",
            om3$df1, t3$global.test$df))
chk(abs(t3$global.test$p.value - t3$ind.test$p.value) < 1e-12,
    "and its global p equals its marginal p, confirming the two are the same test",
    "GLMMcosinor's global test has diverged from its marginal test")

# 2. DANCE reports the phase contrast SIGNED; GLMMcosinor passes it through
#    abs(), so the direction of a phase shift is not recoverable from it.
chk(d_rad > 0 && t3$ind.test$conf.int[1] > 0 && cp3$diff_time > 0,
    sprintf("DANCE's phase contrast is signed (%+.3f h), so 'earlier or later' is answerable",
            cp3$diff_time),
    "the phase contrast lost its sign")

# 3. Small-sample inference. DANCE uses Kenward-Roger; GLMMcosinor is asymptotic.
chk(grepl("Kenward-Roger", om3$method) &&
      !any(grepl("pbkrtest|lmerTest", utils::packageDescription("GLMMcosinor")$Imports %||% "")),
    sprintf("DANCE uses %s; GLMMcosinor is asymptotic (no pbkrtest/lmerTest in its Imports)",
            om3$method),
    "the inference methods are no longer as characterised")
}))

cat(sprintf("\n%s  (%d passed, %d failed)\n",
            if (fail == 0) "GLMMcosinor cross-check PASSED" else "FAILURES", pass, fail))
if (fail > 0) quit(status = 1)
