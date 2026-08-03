# =============================================================================
# test_cmle.R  --  checks for the conditional maximum likelihood estimator
#                  (winstepper_cmle.R)
#
# Run from inside the winstepper folder:   Rscript test_cmle.R
# No Shiny, no network, no API. Pure R.
#
# CML is not a small change to JMLE, it is a different estimator, so these
# checks are the evidence that it is right - not a formality. Tests 1-3 are
# exact: they compare against brute-force enumeration of every response
# pattern, against a closed-form solution, and against a numerically
# differentiated likelihood. Tests 4-6 are statistical recovery checks.
# =============================================================================

source("rasch_engine.R")
source("winsteps_plots.R")
source("winstepper_extras.R")   # the app runs with these overrides in place
source("winstepper_cmle.R")

ok <- 0L; bad <- 0L
check <- function(label, cond, detail = "") {
  if (isTRUE(cond)) { ok <<- ok + 1L; cat(sprintf("  PASS  %s\n", label)) }
  else { bad <<- bad + 1L; cat(sprintf("  FAIL  %s %s\n", label, detail)) }
}

#' Simulate partial-credit responses: one theta per person, its own tau per item.
sim_pcm <- function(n, delta, taus, seed = 1, sd_theta = 1) {
  set.seed(seed)
  th <- stats::rnorm(n, 0, sd_theta)
  I <- length(delta)
  X <- matrix(NA_integer_, n, I,
              dimnames = list(sprintf("P%04d", seq_len(n)), sprintf("I%02d", seq_len(I))))
  for (i in seq_len(I)) {
    tau <- taus[[i]]; m <- length(tau); ct <- c(0, cumsum(tau))
    psi <- matrix(vapply(0:m, function(k) k * (th - delta[i]) - ct[k + 1L], numeric(n)),
                  nrow = n)
    P <- exp(psi - apply(psi, 1L, max)); P <- P / rowSums(P)
    X[, i] <- rowSums(t(apply(P, 1L, cumsum)) < stats::runif(n))
  }
  X
}

#' Dichotomous / rating-scale special case: every item shares one tau.
sim_rasch <- function(n, delta, tau = 0, seed = 1, sd_theta = 1)
  sim_pcm(n, delta, rep(list(tau), length(delta)), seed, sd_theta)

# ---------------------------------------------------------------------------
cat("\n1. EXACT: elementary symmetric functions vs brute-force enumeration\n")
# esf_log() builds gamma_r by convolving items in one at a time. The reference
# sums exp(-sum b) over all prod(m_i + 1) response patterns, grouped by total
# score. Same quantity, completely different route.
loge <- list(-c(0, 0.30, 1.10), -c(0, -0.20, 0.50),
             -c(0, 0.80, 1.40), -c(0, 0.10, -0.30))
a <- esf_log(loge); b <- esf_log_brute(loge)
check("4 items x 3 categories agree to 1e-10", max(abs(a - b)) < 1e-10,
      sprintf("(max |diff| = %.3g)", max(abs(a - b))))
check("gamma_0 = 1 exactly", isTRUE(all.equal(a[1L], 0)))
check("length is sum(m_i) + 1", length(a) == 9L, sprintf("(got %d)", length(a)))

# Unequal category counts, including a dichotomy, and widely spread parameters
# so that the raw gammas would overflow without the log domain.
loge2 <- list(-c(0, 4.0), -c(0, -3.5, 2.0, 9.0), -c(0, 1.0, -6.0),
              -c(0, 0.5, 0.5, 0.5, 12.0))
a2 <- esf_log(loge2); b2 <- esf_log_brute(loge2)
check("mixed category counts and extreme parameters agree to 1e-10",
      max(abs(a2 - b2)) < 1e-10, sprintf("(max |diff| = %.3g)", max(abs(a2 - b2))))
check("no overflow: all log-gammas finite", all(is.finite(a2)))

# ---------------------------------------------------------------------------
cat("\n2. EXACT: two dichotomous items have a closed-form CML solution\n")
# Conditioning on r = 1, P(1,0) = exp(-d1)/(exp(-d1)+exp(-d2)), so the CML
# estimate satisfies  d1 - d2 = log(n01 / n10)  exactly.
n00 <- 20L; n10 <- 30L; n01 <- 15L; n11 <- 25L
X2 <- rbind(matrix(c(0, 0), n00, 2, byrow = TRUE),
            matrix(c(1, 0), n10, 2, byrow = TRUE),
            matrix(c(0, 1), n01, 2, byrow = TRUE),
            matrix(c(1, 1), n11, 2, byrow = TRUE))
colnames(X2) <- c("A", "B"); rownames(X2) <- sprintf("P%03d", seq_len(nrow(X2)))
f2 <- rasch_cmle(rasch_prep(X2, recode = "none"), start = "prox")
want <- log(n01 / n10)
got  <- f2$delta[1L] - f2$delta[2L]
check("delta_A - delta_B matches log(n01/n10) to 1e-6", abs(got - want) < 1e-6,
      sprintf("(got %.8f, expected %.8f)", got, want))
check("items are centred on zero", abs(sum(f2$delta)) < 1e-8)

# ---------------------------------------------------------------------------
cat("\n3. EXACT: analytic gradient vs numerical gradient of the brute-force\n")
cat("   conditional log-likelihood\n")
# The conditional log-likelihood computed by enumerating every pattern shares
# no code with .cml_negll()'s gradient, so this pins the gradient formula
# (d logL/d b_ik = -(N_ik - E[N_ik]) chained through delta and tau).
Xs <- sim_rasch(300, delta = c(-1.0, -0.3, 0.2, 1.1), tau = c(-0.6, 0.6), seed = 7)
ps <- rasch_prep(Xs, recode = "none")
Ig <- ncol(Xs); glg <- "R1"; mccg <- ps$max_cat[glg]; grpg <- ps$groups

brute_cll <- function(par) {
  pr <- .cml_unpack(par, Ig, glg, mccg)
  lg <- .cml_loge(pr$delta, pr$tau, grpg, mccg)
  G  <- esf_log_brute(lg)
  r  <- rowSums(Xs)
  sum(vapply(seq_len(nrow(Xs)), function(n)
    sum(vapply(seq_len(Ig), function(i) lg[[i]][[Xs[n, i] + 1L]], numeric(1))), numeric(1))) -
    sum(G[r + 1L])
}
stats_ <- list(N_ik = .cml_counts(Xs, mccg, grpg), S_i = colSums(Xs),
               C_gj = list(R1 = vapply(1:mccg[[1L]], function(k) sum(Xs >= k), numeric(1))),
               pats = .cml_patterns(Xs, mccg, grpg))
ana <- function(par) -attr(.cml_negll(par, Ig, glg, mccg, grpg, stats_$pats,
                                      stats_$N_ik, stats_$S_i, stats_$C_gj, TRUE), "gradient")
set.seed(11)
p_test <- c(stats::rnorm(Ig - 1L, 0, 0.5), stats::rnorm(mccg[[1L]] - 1L, 0, 0.4))
num <- vapply(seq_along(p_test), function(j) {
  h <- 1e-5; e <- numeric(length(p_test)); e[j] <- h
  (brute_cll(p_test + e) - brute_cll(p_test - e)) / (2 * h)
}, numeric(1))
check("analytic gradient matches the numerical one to 1e-5",
      max(abs(ana(p_test) - num)) < 1e-5,
      sprintf("(max |diff| = %.3g)", max(abs(ana(p_test) - num))))

cat("\n   ...and the CML equations hold at the solution\n")
f3 <- rasch_cmle(ps, conv = 1e-5, start = "prox")
check("estimation converged", isTRUE(f3$converged))
check("observed = conditionally expected sufficient statistic to 1e-5",
      f3$cml$max_residual < 1e-5, sprintf("(got %.3g)", f3$cml$max_residual))
p_hat <- .cml_pack(f3$delta, f3$tau, glg)
check("brute-force numerical gradient is zero at the solution too",
      max(abs(vapply(seq_along(p_hat), function(j) {
        h <- 1e-5; e <- numeric(length(p_hat)); e[j] <- h
        (brute_cll(p_hat + e) - brute_cll(p_hat - e)) / (2 * h)
      }, numeric(1)))) < 1e-3)

# ---------------------------------------------------------------------------
cat("\n4. Parameter recovery from simulated data\n")
true_d <- seq(-1.6, 1.6, length.out = 10)
Xd <- sim_rasch(3000, delta = true_d, seed = 21)
fd <- rasch_cmle(rasch_prep(Xd, recode = "none"))
rd <- stats::cor(fd$delta, true_d); rmse_d <- sqrt(mean((fd$delta - true_d)^2))
cat(sprintf("     dichotomous: r = %.4f, RMSE = %.4f\n", rd, rmse_d))
check("dichotomous item difficulties recovered (r > .99)", rd > 0.99)
check("dichotomous RMSE below 0.10", rmse_d < 0.10, sprintf("(got %.4f)", rmse_d))

true_t <- c(-1.2, -0.1, 1.3)
Xr <- sim_rasch(2500, delta = seq(-1.2, 1.2, length.out = 8), tau = true_t, seed = 31)
fr <- rasch_cmle(rasch_prep(Xr, recode = "none"))
rt <- max(abs(fr$tau$R1 - (true_t - mean(true_t))))
cat(sprintf("     RSM: item r = %.4f, max |tau error| = %.4f\n",
            stats::cor(fr$delta, seq(-1.2, 1.2, length.out = 8)), rt))
check("RSM item difficulties recovered (r > .99)",
      stats::cor(fr$delta, seq(-1.2, 1.2, length.out = 8)) > 0.99)
check("RSM Andrich thresholds recovered to within 0.12", rt < 0.12,
      sprintf("(got %.4f)", rt))
check("RSM standard errors are finite and positive",
      all(is.finite(fr$se_delta)) && all(fr$se_delta > 0) &&
        all(is.finite(fr$se_tau$R1)) && all(fr$se_tau$R1 > 0))

# ---------------------------------------------------------------------------
cat("\n5. PCM: each item its own threshold set\n")
true_dp <- c(-0.8, 0.0, 0.3, 0.9)
Xp <- sim_pcm(2000, delta = true_dp,
              taus = list(c(-1.0, 1.0), c(-1.0, 1.0), c(-0.2, 0.2), c(-0.2, 0.2)),
              seed = 41)
pp <- rasch_prep(Xp, groups = colnames(Xp), recode = "none")
fp <- rasch_cmle(pp)
check("PCM returns one threshold set per item", length(fp$tau) == 4L)
check("PCM converged", isTRUE(fp$converged),
      sprintf("(max residual %.3g)", fp$cml$max_residual))
check("PCM item ordering is recovered",
      stats::cor(fp$delta, true_dp) > 0.95,
      sprintf("(r = %.3f)", stats::cor(fp$delta, true_dp)))
check("PCM thresholds separate the two designs",
      diff(fp$tau[[1L]]) > diff(fp$tau[[3L]]),
      sprintf("(%.2f vs %.2f)", diff(fp$tau[[1L]]), diff(fp$tau[[3L]])))

# ---------------------------------------------------------------------------
cat("\n6. CML vs JMLE: the known JMLE spread inflation\n")
# With a fixed test length L, JMLE item estimates are inflated by roughly
# L/(L-1) - the bias CML does not have. Both estimators run on the same data,
# so the ratio is far more stable than either SD on its own.
L <- 5L
Xb <- sim_rasch(4000, delta = seq(-1.5, 1.5, length.out = L), seed = 51)
pb <- rasch_prep(Xb, recode = "none")
fj <- rasch_jmle(pb, maxit = 500, conv = 1e-5, rconv = 1e-4)
fc <- rasch_cmle(pb)
ratio <- stats::sd(fj$delta) / stats::sd(fc$delta)
cat(sprintf("     L = %d, expected ~%.3f, observed %.3f\n", L, L / (L - 1), ratio))
check("JMLE spread exceeds CML spread", ratio > 1.0, sprintf("(got %.3f)", ratio))
check("ratio is near L/(L-1)", abs(ratio - L / (L - 1)) < 0.15,
      sprintf("(got %.3f, expected %.3f)", ratio, L / (L - 1)))
check("the two estimators order the items identically",
      stats::cor(fj$delta, fc$delta) > 0.999,
      sprintf("(r = %.5f)", stats::cor(fj$delta, fc$delta)))

# ---------------------------------------------------------------------------
cat("\n7. Structure: a CML fit is a drop-in raschfit\n")
fdemo <- rasch_cmle(rasch_prep(as.matrix(
  demo_data(n = 400, k = 8, m = 3)[, paste0("Item", sprintf("%02d", 1:8))]),
  recode = "shift"))
check("class is raschfit", inherits(fdemo, "raschfit"))
for (fld in c("X", "mask", "groups", "max_cat", "cap", "theta", "se_theta",
              "delta", "se_delta", "tau", "se_tau", "E", "W", "C",
              "keep_p", "keep_i", "iterations", "converged", "change_history",
              "settings", "person_id", "item_id"))
  check(sprintf("field '%s' present", fld), !is.null(fdemo[[fld]]))
check("settings record the method", identical(fdemo$settings$method, "CMLE"))
check("every non-extreme person has a measure",
      all(is.finite(fdemo$theta[fdemo$keep_p])))
check("downstream tables run on a CML fit", {
  r <- try({
    item_table(fdemo); person_table(fdemo); summary_table(fdemo)
    category_table(fdemo); score_table(fdemo); pca_residuals(fdemo)
  }, silent = TRUE)
  !inherits(r, "try-error")
})

# ---------------------------------------------------------------------------
cat("\n8. Missing data is conditioned within each pattern\n")
Xm <- sim_rasch(1500, delta = seq(-1.2, 1.2, length.out = 6), seed = 61)
set.seed(62)
Xm[cbind(sample(nrow(Xm), 400, TRUE), sample(ncol(Xm), 400, TRUE))] <- NA
fm <- rasch_cmle(rasch_prep(Xm, recode = "none"))
check("more than one missing-data pattern was found", fm$cml$n_patterns > 1L,
      sprintf("(got %d)", fm$cml$n_patterns))
check("converged with missing data", isTRUE(fm$converged),
      sprintf("(max residual %.3g)", fm$cml$max_residual))
check("difficulties still recovered with ~4%% missing",
      stats::cor(fm$delta, seq(-1.2, 1.2, length.out = 6)) > 0.99)

cat(sprintf("\n%d passed, %d failed\n", ok, bad))
if (bad > 0) quit(status = 1)
