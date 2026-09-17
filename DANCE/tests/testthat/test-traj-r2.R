# ==============================================================================
# tests/testthat/test-traj-r2.R
#
# dance_traj_r2() is an in-house implementation of the Nakagawa & Schielzeth
# (2013) marginal / conditional R-squared with Johnson's (2014) random-slope
# extension. It is checked here against a hand computation from the fitted
# model's own components -- and, when the `performance` package happens to be
# installed, against performance::r2_nakagawa().
# ==============================================================================
.libPaths(c("~/Rlib", .libPaths()))
app_dir <- if (dir.exists("server")) "." else if (dir.exists("../../server")) "../.." else "DANCE"
`%||%` <- function(a, b) if (is.null(a)) b else a
e <- new.env()
for (f in c("08_helpers_cosinor", "08c_helpers_circstat", "08d_helpers_traj",
            "08e_helpers_trajfit", "08f_helpers_trajinf", "08g_helpers_trajband"))
  sys.source(file.path(app_dir, "server", paste0(f, ".R")), envir = e)

mk_long <- function(n = 20, seed = 8) {
  set.seed(seed)
  tp <- seq(0, 22, length.out = 12)
  d <- do.call(rbind, lapply(seq_len(n), function(i) {
    u <- rnorm(1, 0, 3); a <- 7 + rnorm(1, 0, 1.5)
    data.frame(subject = sprintf("S%02d", i), t = tp,
               y = 30 + u + a * cos(2 * pi * (tp - 4) / 24) + rnorm(length(tp), 0, 1.5))
  }))
  d$c1 <- cos(2 * pi * d$t / 24); d$s1 <- sin(2 * pi * d$t / 24)
  d$subject <- factor(d$subject)
  d
}

test_that("random intercept only: Johnson's extension reduces to Nakagawa's formula", {
  skip_if_not_installed("lme4")
  d <- mk_long()
  m <- lme4::lmer(y ~ c1 + s1 + (1 | subject), data = d, REML = TRUE)
  r2 <- e$dance_traj_r2(list(ok = TRUE, model = m))
  expect_false(is.null(r2))
  var_f <- var(as.numeric(predict(m, re.form = NA)))
  var_u <- as.numeric(lme4::VarCorr(m)$subject[1, 1])
  var_e <- sigma(m)^2
  expect_equal(r2$var_fixed, var_f, tolerance = 1e-10)
  expect_equal(r2$var_random, var_u, tolerance = 1e-10)
  expect_equal(r2$var_resid, var_e, tolerance = 1e-10)
  expect_equal(r2$marginal, var_f / (var_f + var_u + var_e), tolerance = 1e-10)
  expect_equal(r2$conditional, (var_f + var_u) / (var_f + var_u + var_e), tolerance = 1e-10)
  expect_true(r2$marginal > 0 && r2$marginal < r2$conditional && r2$conditional < 1)
})

test_that("random slopes: the random-effect variance is the mean of z_i' Sigma z_i", {
  skip_if_not_installed("lme4")
  d <- mk_long()
  m <- lme4::lmer(y ~ c1 + s1 + (1 + c1 + s1 | subject), data = d, REML = TRUE)
  r2 <- e$dance_traj_r2(list(ok = TRUE, model = m))
  expect_false(is.null(r2))
  # hand computation straight from the definition, without Z: every observation's
  # random-effect design row is (1, c1, s1) and Sigma is the subject covariance
  S <- as.matrix(lme4::VarCorr(m)$subject)[c("(Intercept)", "c1", "s1"), c("(Intercept)", "c1", "s1")]
  Zr <- cbind(1, d$c1, d$s1)
  var_r <- mean(rowSums((Zr %*% S) * Zr))
  expect_equal(r2$var_random, var_r, tolerance = 1e-8)
  var_f <- var(as.numeric(predict(m, re.form = NA))); var_e <- sigma(m)^2
  expect_equal(r2$marginal, var_f / (var_f + var_r + var_e), tolerance = 1e-8)
  expect_equal(r2$conditional, (var_f + var_r) / (var_f + var_r + var_e), tolerance = 1e-8)
  # and the diagonal-only shortcut would be WRONG here: it ignores that c1 and
  # s1 are not constant 1, so it is not what the function returns
  expect_false(isTRUE(all.equal(r2$var_random, sum(diag(S)))))
})

test_that("agrees with performance::r2_nakagawa when that package is available", {
  skip_if_not_installed("lme4"); skip_if_not_installed("performance")
  d <- mk_long()
  # the random-intercept model, where performance's own random-effect variance
  # is always available; on a random-slope fit it can decline to compute one
  # ("Random effect variances not available") and then returns a marginal R2
  # that ignores the random part, which is not the quantity under test
  m <- lme4::lmer(y ~ c1 + s1 + (1 | subject), data = d, REML = TRUE)
  r2 <- e$dance_traj_r2(list(ok = TRUE, model = m))
  ref <- suppressWarnings(performance::r2_nakagawa(m))
  skip_if(is.na(as.numeric(ref$R2_conditional)), "performance could not compute the conditional R2")
  expect_equal(unname(r2$marginal), unname(as.numeric(ref$R2_marginal)), tolerance = 1e-3)
  expect_equal(unname(r2$conditional), unname(as.numeric(ref$R2_conditional)), tolerance = 1e-3)
  m2 <- lme4::lmer(y ~ c1 + s1 + (1 + c1 + s1 | subject), data = d, REML = TRUE)
  ref2 <- suppressWarnings(tryCatch(performance::r2_nakagawa(m2), error = function(e) NULL))
  if (!is.null(ref2) && is.finite(as.numeric(ref2$R2_conditional))) {
    r22 <- e$dance_traj_r2(list(ok = TRUE, model = m2))
    expect_equal(unname(r22$marginal), unname(as.numeric(ref2$R2_marginal)), tolerance = 1e-3)
    expect_equal(unname(r22$conditional), unname(as.numeric(ref2$R2_conditional)), tolerance = 1e-3)
  }
})

test_that("a non-merMod fit returns NULL rather than a number", {
  expect_null(e$dance_traj_r2(list(ok = TRUE, model = lm(y ~ c1, data = mk_long()))))
  expect_null(e$dance_traj_r2(list(ok = FALSE)))
})
