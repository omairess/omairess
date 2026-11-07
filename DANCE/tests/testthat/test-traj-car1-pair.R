# ==============================================================================
# tests/testthat/test-traj-car1-pair.R -- independent against CAR(1) residuals
#
# The two residual options of the mixed cosinor must be one pair of models that
# differ in the residual covariance ONLY:
#   * the same observations, fixed-effect design and random structure;
#   * phi = 0 reproduces the independent likelihood at identical parameters;
#   * both REML likelihoods are correct (checked against a dense computation);
#   * the AR(1) optimum is at least the independent one, from either start;
#   * phi is counted as one parameter more (and an estimated tau as one);
#   * the likelihood ratio, AIC and BIC come from the same pair of fits;
#   * the AR(1) fit never moves to another random structure: it refuses, and
#     an estimate bounded by non-convergence is reported as a bound.
# ==============================================================================
.libPaths(c("~/Rlib", .libPaths()))
`%||%` <- function(a, b) if (is.null(a)) b else a
app_dir <- if (dir.exists("server")) "." else if (dir.exists("../../server")) "../.." else "DANCE"
for (f in c("08d_helpers_traj.R", "08e_helpers_trajfit.R", "08f_helpers_trajinf.R",
            "08g_helpers_trajband.R", "08i_helpers_trajsat.R"))
  eval(parse(file.path(app_dir, "server", f), encoding = "UTF-8"), envir = environment())
skip_if_not_installed("lme4")

set.seed(17)
tt <- c(0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 31)
n <- 20; grp <- rep(c("A", "B"), each = 10)
Lc <- t(chol(0.6^abs(outer(tt, tt, "-"))))
Y <- t(vapply(seq_len(n), function(i)
  20 + (grp[i] == "B") + rnorm(1, 0, 1.5) + (4 + rnorm(1, 0, 0.8)) * cos(2 * pi * tt / 24) +
    (2 + rnorm(1, 0, 0.8)) * sin(2 * pi * tt / 24) + 0.05 * tt + as.numeric(Lc %*% rnorm(length(tt))),
  numeric(length(tt))))
Y[3, c(2, 9)] <- NA; Y[14, 12] <- NA
d <- dance_traj_long(Y, tt, sprintf("S%02d", seq_len(n)), list(Group = grp))
roles <- c(Group = "between")
run <- function(res, ...) dance_traj_run(d, period = 24, n_harmonics = 1, roles = roles, t0 = 0,
                                         residuals = res, ...)

# the REML (or ML) log-likelihood written out densely, in lme4's convention
dense_ll <- function(fit, REML = TRUE) {
  m <- fit$model; sp <- fit$spec; dd <- sp$data; phi <- fit$phi %||% 0
  lf <- lme4::lFormula(stats::as.formula(fit$formula), data = dd, REML = REML)
  X <- lf$X; y <- lf$fr[[1]]; Lam <- lf$reTrms$Lambdat
  Lam@x <- lme4::getME(m, "theta")[lf$reTrms$Lind]
  cv <- as.character(dd[[sp$curve_var %||% "subject"]])
  R <- if (phi > 0) outer(seq_len(nrow(dd)), seq_len(nrow(dd)), function(i, j)
    ifelse(cv[i] == cv[j], phi^abs(dd$t[i] - dd$t[j]), 0)) else diag(nrow(dd))
  ZL <- as.matrix(Matrix::t(lf$reTrms$Zt) %*% Matrix::t(Lam))
  V <- stats::sigma(m)^2 * (ZL %*% t(ZL) + R); Vi <- solve(V)
  A <- t(X) %*% Vi %*% X; r <- y - X %*% solve(A, t(X) %*% Vi %*% y)
  q <- as.numeric(t(r) %*% Vi %*% r); ldV <- as.numeric(determinant(V)$modulus)
  if (!REML) return(-0.5 * (length(y) * log(2 * pi) + ldV + q))
  -0.5 * ((length(y) - ncol(X)) * log(2 * pi) + ldV + as.numeric(determinant(A)$modulus) + q)
}

ind <- run("independent", trend = "linear")
car <- run("car1", trend = "linear")

test_that("both options use the same observations, fixed design and random structure", {
  expect_true(ind$ok && car$ok)
  expect_identical(stats::nobs(ind$model), stats::nobs(car$model))
  expect_identical(ind$spec$data[c("subject", "t", "y")], car$spec$data[c("subject", "t", "y")])
  expect_identical(ind$spec$fixed_formula, car$spec$fixed_formula)
  Xc <- lme4::lFormula(stats::as.formula(car$formula), data = car$spec$data)$X
  expect_equal(c(as.matrix(lme4::getME(ind$model, "X"))), c(as.matrix(Xc)))
  W <- dance_car1_whitener(car$spec$data$t, car$spec$data$subject, car$phi)
  expect_equal(c(as.matrix(W %*% Xc)), c(as.matrix(lme4::getME(car$model, "X"))))
  expect_identical(ind$re_formula, car$re_formula)
  expect_identical(ind$re_rung, car$re_rung)
})

test_that("phi = 0 reproduces the independent likelihood at identical parameters", {
  lf <- lme4::lFormula(stats::as.formula(ind$formula), data = ind$spec$data, REML = TRUE)
  dv_i <- do.call(lme4::mkLmerDevfun, lf)
  th <- lme4::getME(ind$model, "theta")
  for (ph in c(0, 1e-8)) {
    lw <- lf; Wp <- dance_car1_whitener(ind$spec$data$t, ind$spec$data$subject, ph)
    att <- attributes(lw$X); Xw <- as.matrix(Wp %*% lw$X); attributes(Xw) <- att; lw$X <- Xw
    lw$fr[[1]] <- as.numeric(Wp %*% lw$fr[[1]])
    lw$reTrms$Zt <- methods::as(lw$reTrms$Zt %*% Matrix::t(Wp), "CsparseMatrix")
    dv_w <- do.call(lme4::mkLmerDevfun, lw)
    expect_equal(dv_w(th) - 2 * attr(Wp, "logdet"), dv_i(th), tolerance = 1e-8)
  }
  ctrl <- lme4::lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))
  m0 <- dance_traj_lmer_car1(stats::as.formula(ind$formula), ind$spec$data, 0, "subject", TRUE, ctrl, light = TRUE)
  expect_equal(as.numeric(stats::logLik(m0)), as.numeric(stats::logLik(ind$model)), tolerance = 1e-8)
})

test_that("both REML likelihoods are the data's, as a dense computation gives them", {
  expect_true(dance_traj_loglik(ind)$REML && dance_traj_loglik(car)$REML)
  expect_equal(dance_traj_loglik(ind)$logLik, dense_ll(ind), tolerance = 1e-7)
  expect_equal(dance_traj_loglik(car)$logLik, dense_ll(car), tolerance = 1e-7)
})

test_that("the AR(1) optimum is at least the independent one, from either start", {
  expect_gte(dance_traj_loglik(car)$logLik, dance_traj_loglik(ind)$logLik)
  expect_false(isTRUE(car$residual_ar$optimum_below_independent))
  # started at the independent estimates with phi near zero: the same optimum
  ctrl <- lme4::lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))
  for (ph in c(1e-3, car$phi)) {
    lw <- lme4::lFormula(stats::as.formula(car$formula), data = car$spec$data, REML = TRUE, control = ctrl)
    Wp <- dance_car1_whitener(car$spec$data$t, car$spec$data$subject, ph)
    att <- attributes(lw$X); Xw <- as.matrix(Wp %*% lw$X); attributes(Xw) <- att; lw$X <- Xw
    lw$fr[[1]] <- as.numeric(Wp %*% lw$fr[[1]])
    lw$reTrms$Zt <- methods::as(lw$reTrms$Zt %*% Matrix::t(Wp), "CsparseMatrix")
    dv <- do.call(lme4::mkLmerDevfun, c(lw, list(control = ctrl)))
    o1 <- lme4::optimizeLmer(dv, optimizer = "bobyqa", control = list(maxfun = 2e5))
    th <- lme4::getME(ind$model, "theta")
    o2 <- lme4::optimizeLmer(dv, optimizer = "bobyqa", start = list(theta = th), control = list(maxfun = 2e5))
    expect_equal(o1$fval, o2$fval, tolerance = 1e-6)
    # and no point of the profile lies above the reported maximum
    expect_lte(max(car$residual_ar$profile$logLik), dance_traj_loglik(car)$logLik + 1e-6)
  }
})

test_that("phi is counted as one parameter more", {
  li <- dance_traj_loglik(ind); lc <- dance_traj_loglik(car)
  p <- ncol(lme4::getME(ind$model, "X")); th <- length(lme4::getME(ind$model, "theta"))
  expect_identical(li$df, as.integer(p + th + 1L))
  expect_identical(lc$df, li$df + 1L)
})

test_that("the likelihood ratio, AIC and BIC come from the same pair of fitted models", {
  ct <- car$residual_ar$comparison
  li <- dance_traj_loglik(ind); lc <- dance_traj_loglik(car)
  expect_identical(ct$likelihood, c("REML", "REML"))
  expect_equal(ct$logLik, c(li$logLik, lc$logLik))
  expect_equal(ct$AIC, c(li$AIC, lc$AIC)); expect_equal(ct$BIC, c(li$BIC, lc$BIC))
  expect_identical(ct$random_effects[1], ct$random_effects[2])
  expect_equal(car$residual_ar$lrt$statistic, 2 * (lc$logLik - li$logLik))
  expect_equal(car$residual_ar$lrt$p, 0.5 * stats::pchisq(2 * (lc$logLik - li$logLik), 1, lower.tail = FALSE))
})

test_that("with tau estimated, the comparison is by ML from one refitted pair, tau counted", {
  run_t <- function(res) run(res, trend = "exp_sat", tau_mode = "estimate", tau_value = 14,
                             tau_lo = 4, tau_hi = 40)
  ci <- run_t("independent"); cc <- run_t("car1")
  expect_true(ci$ok && cc$ok)
  expect_identical(ci$re_formula, cc$re_formula)
  ct <- cc$residual_ar$comparison
  expect_identical(ct$likelihood, c("ML", "ML"))
  expect_identical(ct$df[2], ct$df[1] + 1L)
  # tau is a parameter of both: fixed effects + variance components + sigma + tau (+ phi)
  p <- ncol(lme4::getME(ci$model, "X")); th <- length(lme4::getME(ci$model, "theta"))
  expect_identical(ct$df[1], as.integer(p + th + 2L))
  expect_equal(cc$residual_ar$lrt$statistic, max(0, 2 * diff(ct$logLik)))
  expect_gte(ct$logLik[2], ct$logLik[1])
  # the fitted model's own criteria count the estimated tau too
  expect_identical(dance_traj_loglik(ci)$df, as.integer(p + th + 2L))
})

test_that("the AR(1) fit never moves to another random structure", {
  orig <- dance_traj_fit_one
  on.exit(assign("dance_traj_fit_one", orig, envir = parent.frame()), add = TRUE)
  # the final refit at phi-hat fails on the independent model's structure
  dance_traj_fit_one <<- function(spec, re_formula, engine = "lmer", REML = TRUE, cor_struct = NULL,
                                  fixed_formula = NULL, light = FALSE, phi = NULL) {
    if (!isTRUE(light) && isTRUE(phi > 0.3) && identical(re_formula, spec$re_ladder[[1]]$formula)) return(NULL)
    orig(spec, re_formula, engine, REML, cor_struct, fixed_formula, light, phi)
  }
  r <- run("car1")
  expect_false(r$ok)
  expect_match(r$message, "random structure the independent-residual model uses", fixed = TRUE)
  # every fit above phi = 0.3 fails: the estimate is a bound, and says so
  dance_traj_fit_one <<- function(spec, re_formula, engine = "lmer", REML = TRUE, cor_struct = NULL,
                                  fixed_formula = NULL, light = FALSE, phi = NULL) {
    if (isTRUE(phi > 0.3) && identical(re_formula, spec$re_ladder[[1]]$formula)) return(NULL)
    orig(spec, re_formula, engine, REML, cor_struct, fixed_formula, light, phi)
  }
  r <- run("car1")
  expect_true(r$ok)
  expect_identical(r$re_rung, 1L)
  expect_true(r$residual_ar$bounded_by_failure)
  expect_false(r$residual_ar$at_upper)
  expect_match(dance_traj_residual_statement(r), "a bound set by the optimiser", fixed = TRUE)
  dance_traj_fit_one <<- orig
})
