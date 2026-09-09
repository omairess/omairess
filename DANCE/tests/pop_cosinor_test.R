# ==============================================================================
# tests/pop_cosinor_test.R — the population-mean cosinor (Bingham et al. 1982)
#
# The formulas came out of a reference implementation's documentation, and a
# transcribed formula is exactly the kind of thing that goes wrong silently:
# a sign on a cross term or a wrong denominator degrees of freedom shifts the
# Type I error without producing anything that looks broken. So the test is
# CALIBRATION, not a comparison against the transcription.
#
# Under each null in turn, the rejection rate must sit at the nominal level. The
# amplitude and acrophase tests are the informative ones, because those are the
# two whose denominators depend on the pooled mean direction -- the part most
# easily got wrong.
#
# Run with:   Rscript tests/pop_cosinor_test.R      (from the DANCE directory)
# ==============================================================================
.libPaths(c("~/Rlib", .libPaths()))
ok_n <- 0L; bad <- 0L
chk <- function(cond, good, bad_msg) {
  if (isTRUE(cond)) { cat("ok   ", good, "\n"); ok_n <<- ok_n + 1L }
  else { cat("FAIL:", bad_msg, "\n"); bad <<- bad + 1L }
}
app_dir <- if (dir.exists("server")) "." else "DANCE"
e <- new.env(parent = globalenv())
source(file.path(app_dir, "server/08b_helpers_popcosinor.R"), local = e)

# One simulated study: m groups of n participants, each contributing a (cos, sin)
# coefficient pair and a MESOR drawn around that group's true values.
sim <- function(m = 3, n = 15, amp = c(10, 10, 10), acro = c(2, 2, 2),
                mes = c(50, 50, 50), sd_coef = 3, sd_mes = 4) {
  g <- rep(sprintf("g%d", seq_len(m)), each = n)
  b <- unlist(lapply(seq_len(m), function(j) amp[j] * cos(acro[j]) + rnorm(n, 0, sd_coef)))
  s <- unlist(lapply(seq_len(m), function(j) amp[j] * sin(acro[j]) + rnorm(n, 0, sd_coef)))
  M <- unlist(lapply(seq_len(m), function(j) rnorm(n, mes[j], sd_mes)))
  list(beta = b, gamma = s, mesor = M, group = factor(g))
}

rate <- function(param, nsim = 2000, ...) {
  p <- vapply(seq_len(nsim), function(i) {
    set.seed(5000 + i)
    d <- sim(...)
    r <- e$dance_pop_cosinor(d$beta, d$gamma, d$mesor, d$group)
    r$tests$p[r$tests$parameter == param]
  }, numeric(1))
  c(rate = mean(p < .05, na.rm = TRUE),
    ks = suppressWarnings(stats::ks.test(p[is.finite(p)], "punif")$p.value))
}

cat("-- calibration under the null (nominal .05, 2000 simulations) ----------\n")
for (prm in c("MESOR", "Amplitude", "Acrophase")) {
  r <- rate(prm)
  chk(abs(r[["rate"]] - .05) < .02 && r[["ks"]] > .01,
      sprintf("%-9s rejects at %.3f under its null (KS p = %.3f)", prm, r[["rate"]], r[["ks"]]),
      sprintf("%-9s is MIScalibrated: %.3f (KS p = %.3f) -- check the formula and df",
              prm, r[["rate"]], r[["ks"]]))
}

cat("\n-- power: each test must find its own effect and not the others -------\n")
one <- function(...) { set.seed(11); d <- sim(...)
  e$dance_pop_cosinor(d$beta, d$gamma, d$mesor, d$group)$tests }

t_m <- one(mes = c(50, 58, 66))                       # MESOR differs only
chk(t_m$p[1] < .01, sprintf("a MESOR difference is found (p = %.4g)", t_m$p[1]),
    sprintf("a MESOR difference was missed (p = %.4g)", t_m$p[1]))
chk(t_m$p[2] > .05 && t_m$p[3] > .05,
    "a MESOR difference does not leak into the amplitude or acrophase tests",
    sprintf("a MESOR shift contaminated amplitude (p = %.3g) or acrophase (p = %.3g)",
            t_m$p[2], t_m$p[3]))

t_a <- one(amp = c(6, 14, 22))                        # amplitude differs only
chk(t_a$p[2] < .01, sprintf("an amplitude difference is found (p = %.4g)", t_a$p[2]),
    sprintf("an amplitude difference was missed (p = %.4g)", t_a$p[2]))

t_p <- one(acro = c(1.2, 2.2, 3.2))                   # acrophase differs only
chk(t_p$p[3] < .01, sprintf("an acrophase difference is found (p = %.4g)", t_p$p[3]),
    sprintf("an acrophase difference was missed (p = %.4g)", t_p$p[3]))

cat("\n-- Bingham's caution, carried as a value ------------------------------\n")
# Asserting the flag's value on ONE null draw is a flaky test, and was: seed 3
# gave acrophase p = .030, which is a real 5%-rate false positive, not a defect.
# The two things worth asserting are the RELATIONSHIP, which must hold on every
# draw, and the RATE, which must sit at the nominal level.
set.seed(3); d0 <- sim()
r0 <- e$dance_pop_cosinor(d0$beta, d0$gamma, d0$mesor, d0$group)
chk(identical(r0$amplitude_interpretable, !r0$acrophase_differs),
    "the amplitude caution is exactly the negation of the acrophase verdict",
    "the amplitude flag and the acrophase verdict disagree")

set.seed(21); dA <- sim(acro = c(1.2, 2.2, 3.2))
rA <- e$dance_pop_cosinor(dA$beta, dA$gamma, dA$mesor, dA$group)
chk(isTRUE(rA$acrophase_differs) && isFALSE(rA$amplitude_interpretable),
    "with genuinely different acrophases, the amplitude comparison is flagged",
    "a real acrophase difference did not raise Bingham's amplitude caution")

fires <- mean(vapply(seq_len(300), function(i) {
  set.seed(9000 + i); dd <- sim()
  isTRUE(e$dance_pop_cosinor(dd$beta, dd$gamma, dd$mesor, dd$group)$acrophase_differs)
}, logical(1)))
chk(fires < .10,
    sprintf("under the null the caution fires %.1f%% of the time, near its nominal 5%%",
            100 * fires),
    sprintf("the caution fires %.1f%% of the time under the null", 100 * fires))

cat("\n-- refusals and bookkeeping ------------------------------------------\n")
one_g <- e$dance_pop_cosinor(rnorm(10), rnorm(10), rnorm(10), factor(rep("a", 10)))
chk(!isTRUE(one_g$ok), "a single group is refused", "a single group was accepted")
tiny <- e$dance_pop_cosinor(rnorm(3), rnorm(3), rnorm(3), factor(c("a","a","b")))
chk(!isTRUE(tiny$ok) && grepl("at least 2 participants", tiny$message),
    "a group of one is refused, naming the sizes",
    "a group with one participant was accepted, so the pooled covariance is undefined")
chk(all(r0$tests$df1 == r0$n_groups - 1) && all(r0$tests$df2 == r0$n_subjects - r0$n_groups),
    sprintf("degrees of freedom are (m-1, K-m) = (%d, %d)", r0$df1, r0$df2),
    "the degrees of freedom are not (m-1, K-m)")
chk(all(r0$groups$acrophase_time >= 0 & r0$groups$acrophase_time < 24),
    "acrophases are reported in time units inside one period",
    "an acrophase falls outside the period")

cat(sprintf("\n%s  (%d passed, %d failed)\n",
            if (bad == 0) "Population cosinor tests PASSED" else "Population cosinor tests FAILED",
            ok_n, bad))
quit(status = if (bad == 0) 0 else 1)
