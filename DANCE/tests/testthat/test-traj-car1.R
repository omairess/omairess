# ==============================================================================
# tests/testthat/test-traj-car1.R -- autocorrelated residuals, by whitening
#
# The mixed cosinor's continuous-time AR(1) residual model is fitted as an
# ordinary lme4 model on whitened data, so its Kenward-Roger and Satterthwaite
# tests apply unchanged. These pin the identity that makes that legitimate,
# against nlme's own implementation of the same model (lme with corCAR1):
#   * the whitener W satisfies W R W' = I for uneven gaps and missing times;
#   * at a given phi, the whitened fit IS the corCAR1 model: same fixed
#     effects, standard errors, sigma and REML log-likelihood;
#   * the phi profile finds nlme's REML estimate;
#   * phi = 0 is the independent-residual fit exactly;
#   * the tests run on it, with the df the design implies.
# ==============================================================================
.libPaths(c("~/Rlib", .libPaths()))
`%||%` <- function(a, b) if (is.null(a)) b else a
app_dir <- if (dir.exists("server")) "." else if (dir.exists("../../server")) "../.." else "DANCE"
for (f in c("08d_helpers_traj.R", "08e_helpers_trajfit.R", "08f_helpers_trajinf.R",
            "08g_helpers_trajband.R", "08i_helpers_trajsat.R"))
  eval(parse(file.path(app_dir, "server", f), encoding = "UTF-8"), envir = environment())
skip_if_not_installed("nlme"); skip_if_not_installed("lme4"); skip_if_not_installed("pbkrtest")

# two groups x 10 participants, uneven times, a few missing, CAR(1) residuals
set.seed(17)
tt <- c(0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 31)
n <- 20; grp <- rep(c("A", "B"), each = 10)
Lc <- t(chol(0.6^abs(outer(tt, tt, "-"))))
Y <- t(vapply(seq_len(n), function(i)
  20 + (grp[i] == "B") + rnorm(1, 0, 1.5) + (4 + rnorm(1, 0, 0.8)) * cos(2 * pi * tt / 24) +
    (2 + rnorm(1, 0, 0.8)) * sin(2 * pi * tt / 24) + as.numeric(Lc %*% rnorm(length(tt))),
  numeric(length(tt))))
Y[3, c(2, 9)] <- NA; Y[14, 12] <- NA
d <- dance_traj_long(Y, tt, sprintf("S%02d", seq_len(n)), list(Group = grp))
sp <- dance_traj_spec(d, period = 24, n_harmonics = 1, trend = "none", t0 = 0,
                      roles = c(Group = "between"))

test_that("the whitener turns the CAR(1) correlation into the identity", {
  dd <- sp$data; cv <- as.character(dd$curve)
  W <- as.matrix(dance_car1_whitener(dd$t, cv, 0.6))
  R <- outer(seq_len(nrow(dd)), seq_len(nrow(dd)), function(i, j)
    ifelse(cv[i] == cv[j], 0.6^abs(dd$t[i] - dd$t[j]), 0))
  expect_equal(W %*% R %*% t(W), diag(nrow(dd)), tolerance = 1e-10)
  # log|det W| is minus half log|det R|
  expect_equal(attr(dance_car1_whitener(dd$t, cv, 0.6), "logdet"),
               -0.5 * as.numeric(determinant(R, logarithm = TRUE)$modulus), tolerance = 1e-8)
  # phi = 0 is no transformation at all
  expect_equal(as.matrix(dance_car1_whitener(dd$t, cv, 0)), diag(nrow(dd)))
  expect_error(dance_car1_whitener(c(0, 1, 1), c("a", "a", "a"), 0.5), "share a time point")
})

ref <- nlme::lme(y ~ (c1 + s1) * Group, random = ~ 1 + c1 + s1 | subject, data = sp$data,
                 correlation = nlme::corCAR1(form = ~ t | subject), method = "REML",
                 control = nlme::lmeControl(maxIter = 500, msMaxIter = 500, opt = "optim"))
phi_ref <- as.numeric(stats::coef(ref$modelStruct$corStruct, unconstrained = FALSE))

test_that("at a given phi the whitened fit is nlme's corCAR1 model", {
  fw <- dance_traj_fit(sp, REML = TRUE, max_rung = 1L, phi = phi_ref)
  expect_true(fw$ok)
  expect_equal(fw$re_formula, "(1 + c1 + s1 | subject)")
  b <- lme4::fixef(fw$model)
  expect_equal(unname(b[names(nlme::fixef(ref))]), unname(nlme::fixef(ref)), tolerance = 1e-5)
  se_w <- sqrt(diag(as.matrix(stats::vcov(fw$model))))[names(nlme::fixef(ref))]
  expect_equal(unname(se_w), unname(sqrt(diag(stats::vcov(ref)))), tolerance = 1e-3)
  expect_equal(stats::sigma(fw$model), ref$sigma, tolerance = 1e-4)
  expect_equal(dance_traj_loglik(fw)$logLik, as.numeric(stats::logLik(ref)), tolerance = 1e-5)
})

test_that("the phi profile finds nlme's REML estimate, and tests independence", {
  pp <- dance_traj_profile_phi(sp, rung = 1L, REML = TRUE)
  expect_true(pp$ok)
  expect_equal(pp$phi, phi_ref, tolerance = 2e-3)
  expect_true(pp$ci[1] < pp$phi && pp$ci[2] > pp$phi)
  expect_true(pp$lrt$p < 0.001)
  # the likelihood ratio is the one between nlme's CAR(1) and independent fits
  ref0 <- nlme::lme(y ~ (c1 + s1) * Group, random = ~ 1 + c1 + s1 | subject, data = sp$data,
                    method = "REML",
                    control = nlme::lmeControl(maxIter = 500, msMaxIter = 500, opt = "optim"))
  expect_equal(pp$lrt$statistic,
               2 * (as.numeric(stats::logLik(ref)) - as.numeric(stats::logLik(ref0))),
               tolerance = 0.05)
})

test_that("phi = 0 is the independent-residual fit exactly", {
  f0 <- dance_traj_fit(sp, REML = TRUE, max_rung = 1L, phi = 0)
  fi <- dance_traj_fit(sp, REML = TRUE, max_rung = 1L)
  expect_equal(lme4::fixef(f0$model), lme4::fixef(fi$model))
  expect_equal(dance_traj_loglik(f0)$logLik, as.numeric(stats::logLik(fi$model)))
})

test_that("Kenward-Roger and Satterthwaite run on the whitened fit", {
  fw <- dance_traj_fit(sp, REML = TRUE, max_rung = 1L, phi = phi_ref)
  kr <- dance_traj_marginal_test(fw, "Group", "circadian", df_method = "kr")
  st <- dance_traj_marginal_test(fw, "Group", "circadian", df_method = "satterthwaite")
  expect_true(kr$ok && st$ok)
  expect_match(kr$method, "Kenward-Roger"); expect_match(st$method, "Satterthwaite")
  # 20 participants in 2 groups, random rhythm per participant: about 18 df,
  # where nlme's own containment rule gives the within-participant count
  expect_true(kr$df2 > 15 && kr$df2 < 20)
  expect_true(st$df2 > 15 && st$df2 < 20)
  expect_gt(stats::anova(ref, L = c("c1:GroupB" = 1))$denDF, 200)
})

test_that("the run function estimates phi, refits there, and reports it", {
  run <- dance_traj_run(d, period = 24, n_harmonics = 1, trend = "none",
                        roles = c(Group = "between"), t0 = 0, residuals = "car1")
  expect_true(run$ok)
  expect_identical(run$residual_ar$method, "REML")
  expect_equal(run$phi, run$residual_ar$phi)
  expect_match(dance_traj_residual_statement(run), "continuous-time AR(1)", fixed = TRUE)
  ind <- dance_traj_run(d, period = 24, n_harmonics = 1, trend = "none",
                        roles = c(Group = "between"), t0 = 0)
  expect_identical(ind$phi, 0)
  expect_match(dance_traj_residual_statement(ind), "independent", fixed = TRUE)
})

# ---- a within-participant factor: the AR(1) runs within each CURVE ------------
# two sessions per participant, each its own series: the correlation is by
# curve (participant x session), never across a participant's sessions
set.seed(23)
tw <- c(0, 2, 4, 6, 8, 11, 14, 16, 18, 21)
rw <- expand.grid(i = 1:12, ses = c("A", "B"))
uw <- matrix(rnorm(36, 0, c(1.5, 0.8, 0.8)), 12, 3, byrow = TRUE)
Lw <- t(chol(0.5^abs(outer(tw, tw, "-"))))
Yw <- t(vapply(seq_len(nrow(rw)), function(k) {
  i <- rw$i[k]
  20 + (rw$ses[k] == "B") * 0.8 + uw[i, 1] + rnorm(1, 0, 0.5) +
    (4 + uw[i, 2]) * cos(2 * pi * tw / 24) + (2 + uw[i, 3]) * sin(2 * pi * tw / 24) +
    as.numeric(Lw %*% rnorm(length(tw)))
}, numeric(length(tw))))
Yw[3, 4] <- NA; Yw[15, c(2, 9)] <- NA
dw <- dance_traj_long(Yw, tw, sprintf("S%02d", rw$i), list(Session = as.character(rw$ses)))
spw <- dance_traj_spec(dw, period = 24, n_harmonics = 1, trend = "none", t0 = 0,
                       roles = c(Session = "within"))
ctl_w <- nlme::lmeControl(maxIter = 500, msMaxIter = 500, opt = "optim")

test_that("with a within factor the whitened fit is nlme's corCAR1 by curve", {
  expect_true(spw$has_curve_level)
  # rung 2: participant rhythm plus a curve level; rung 4: random effects by
  # participant only, the correlation still by curve (a finer grouping)
  refs <- list(
    `2` = nlme::lme(y ~ (c1 + s1) * Session, random = list(subject = ~ 1 + c1 + s1, curve = ~ 1),
                    data = spw$data, correlation = nlme::corCAR1(form = ~ t | subject/curve),
                    method = "REML", control = ctl_w),
    `4` = nlme::lme(y ~ (c1 + s1) * Session, random = ~ 1 + c1 + s1 | subject,
                    data = spw$data, correlation = nlme::corCAR1(form = ~ t | subject/curve),
                    method = "REML", control = ctl_w))
  for (k in names(refs)) {
    ref_k <- refs[[k]]
    phi_k <- as.numeric(stats::coef(ref_k$modelStruct$corStruct, unconstrained = FALSE))
    fk <- dance_traj_fit(spw, REML = TRUE, min_rung = as.integer(k), max_rung = as.integer(k),
                         phi = phi_k)
    expect_true(fk$ok)
    expect_equal(fk$re_rung, as.integer(k))
    expect_equal(unname(lme4::fixef(fk$model)[names(nlme::fixef(ref_k))]),
                 unname(nlme::fixef(ref_k)), tolerance = 1e-5)
    expect_equal(dance_traj_loglik(fk)$logLik, as.numeric(stats::logLik(ref_k)), tolerance = 1e-5)
    expect_equal(dance_traj_profile_phi(spw, rung = as.integer(k), REML = TRUE)$phi, phi_k,
                 tolerance = 2e-3)
  }
})

test_that("a fit started below the first rung records the rungs it skipped", {
  # the tau and phi profiles fit every point on the structure already chosen
  # (min_rung); this used to fail outright whenever that was not rung 1
  f3 <- dance_traj_fit(spw, REML = TRUE, min_rung = 3L, max_rung = 3L)
  expect_true(f3$ok)
  expect_equal(f3$re_rung, 3L)
  expect_true(isTRUE(f3$attempts[[1]]$skipped) && isTRUE(f3$attempts[[2]]$skipped))
  expect_length(f3$simplified_reason, 2L)
  expect_match(dance_traj_fit_report(f3), "not refitted (structure carried over)", fixed = TRUE)
})

test_that("p-values format from whatever arrives", {
  expect_identical(dance_fmt_p(0.0234), ".023")
  expect_identical(dance_fmt_p(0.0001), "< .001")
  expect_identical(dance_fmt_p(NA), "—")
  expect_identical(dance_fmt_p(NULL), "—")
  expect_identical(dance_fmt_p(list(0.5)), ".500")
  expect_identical(dance_fmt_p(numeric(0)), "—")
})
