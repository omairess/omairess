# ==============================================================================
# tests/testthat/test-traj-sat.R
#
# The saturating trend's time constant: profiled once for all cells in the
# trajectory model, and estimated per cell -- with A and A/tau -- by a nonlinear
# mixed model (server/08i). These pin the per-cell model to nlme's own numbers
# (its anova F, its likelihood-ratio test), the delta method to a numerical
# derivative, the curves to their closed form, and the profile to one random
# structure with a refined maximum inside its own interval.
# ==============================================================================
.libPaths(c("~/Rlib", .libPaths()))
app_dir <- if (dir.exists("server")) "." else if (dir.exists("../../server")) "../.." else "DANCE"
`%||%` <- function(a, b) if (is.null(a)) b else a

e <- new.env()
for (f in c("08_helpers_cosinor", "08c_helpers_circstat", "08d_helpers_traj",
            "08e_helpers_trajfit", "08f_helpers_trajinf", "08g_helpers_trajband",
            "08h_helpers_twostage", "08i_helpers_trajsat"))
  sys.source(file.path(app_dir, "server", paste0(f, ".R")), envir = e)

clock <- c(8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30, 32, 33)
mk_long <- function(taus = c(A = 5, B = 12), A = c(A = 20, B = 20), nper = 14, seed = 3) {
  set.seed(seed)
  tt <- clock - clock[1]
  grp <- rep(names(taus), each = nper)
  Y <- t(sapply(seq_along(grp), function(i) {
    g <- grp[i]
    30 + rnorm(1, 0, 1.5) + (A[[g]] + rnorm(1, 0, 1)) * (1 - exp(-tt / taus[[g]])) +
      5 * cos(2 * pi * (clock - 16) / 24) + rnorm(length(tt), 0, 1)
  }))
  e$dance_traj_long(Y, clock, sprintf("S%02d", seq_along(grp)), list(Group = grp))
}
mk_fit <- function(d, tau = 8) e$dance_traj_fit(e$dance_traj_spec(d, 24, 1, "exp_sat", tau = tau))

skip_heavy <- function() {
  skip_if_not_installed("lme4"); skip_if_not_installed("nlme")
}

test_that("the per-cell model recovers each cell's A and tau, and knows its df", {
  skip_heavy()
  ff <- mk_fit(mk_long()); skip_if_not(isTRUE(ff$ok))
  sf <- e$dance_traj_sat_fit(ff)
  expect_true(isTRUE(sf$ok))
  expect_equal(sf$n_cells, 2L)
  expect_equal(sf$df, 28L - 2L)                 # participants minus between cells
  st <- sf$table
  expect_equal(st$cell, c("A", "B"))
  # planted tau 5 and 12, A 20: inside their intervals and near the truth
  expect_true(st$tau_lo[1] < 5 && 5 < st$tau_hi[1])
  expect_true(st$tau_lo[2] < 12 && 12 < st$tau_hi[2])
  expect_lt(abs(st$A[1] - 20), 3)
  expect_equal(st$init_slope, st$A / st$tau, tolerance = 1e-12)
  # the interval for tau is the back-transformed log-scale interval: geometric
  # symmetry about the estimate
  expect_equal(log(st$tau) - log(st$tau_lo), log(st$tau_hi) - log(st$tau), tolerance = 1e-10)
})

test_that("its tests are nlme's own: the Wald F equals anova(L =), the LRT equals anova()", {
  skip_heavy()
  ff <- mk_fit(mk_long()); skip_if_not(isTRUE(ff$ok))
  sf <- e$dance_traj_sat_fit(ff); skip_if_not(isTRUE(sf$ok))
  fe <- nlme::fixef(sf$fit)
  L <- matrix(0, 1, length(fe), dimnames = list(NULL, names(fe)))
  L[1, "lt.cell_g1"] <- -1; L[1, "lt.cell_g2"] <- 1
  an <- nlme::anova.lme(sf$fit, L = L)
  w <- e$dance_traj_sat_test(sf, ff, "tau")
  expect_equal(w$statistic, an[["F-value"]], tolerance = 1e-8)
  expect_equal(w$df1, 1)
  expect_equal(w$df2, sf$df)                    # the between-participant df, not nlme's
  # the one-tau test is the likelihood ratio nlme itself reports
  lr <- stats::anova(sf$fit_common, sf$fit)
  expect_equal(sf$lrt$statistic, lr$L.Ratio[2], tolerance = 1e-8)
  expect_equal(sf$lrt$df, 1L)
  expect_gt(sf$lrt$statistic, 0)
})

test_that("A/tau goes through the delta method, checked against a numerical derivative", {
  skip_heavy()
  ff <- mk_fit(mk_long()); skip_if_not(isTRUE(ff$ok))
  sf <- e$dance_traj_sat_fit(ff); skip_if_not(isTRUE(sf$ok))
  V <- e$dance_traj_sat_V(sf$fit)
  for (k in 1:2) {
    nm <- c(e$dance_traj_sat_coef(sf, "A", k), e$dance_traj_sat_coef(sf, "lt", k))
    th <- nlme::fixef(sf$fit)[nm]
    g <- function(p) p[1] / exp(p[2]); h <- 1e-6
    J <- c((g(th + c(h, 0)) - g(th - c(h, 0))) / (2 * h), (g(th + c(0, h)) - g(th - c(0, h))) / (2 * h))
    expect_equal(sf$table$init_slope_se[k], sqrt(drop(t(J) %*% V[nm, nm] %*% J)), tolerance = 1e-6)
  }
  # pairwise: tau as a ratio, A/tau as a difference, each with its own interval
  pt <- e$dance_traj_sat_pairwise(sf, "tau", adjust = "none")
  expect_true(pt$ratio)
  expect_equal(pt$table$estimate, sf$table$tau[1] / sf$table$tau[2], tolerance = 1e-10)
  expect_true(pt$table$lo < pt$table$estimate && pt$table$estimate < pt$table$hi)
  ps <- e$dance_traj_sat_pairwise(sf, "init_slope", adjust = "none")
  expect_false(ps$ratio)
  expect_equal(ps$table$estimate, sf$table$init_slope[1] - sf$table$init_slope[2], tolerance = 1e-10)
})

test_that("the per-cell trend curves are A (1 - exp(-t/tau)), zero at the first observation", {
  skip_heavy()
  ff <- mk_fit(mk_long()); skip_if_not(isTRUE(ff$ok))
  sf <- e$dance_traj_sat_fit(ff); skip_if_not(isTRUE(sf$ok))
  cv <- e$dance_traj_sat_curves(sf, n_time = 11)
  for (k in 1:2) {
    d <- cv$table[cv$table$cell == sf$cells[k], ]
    expect_equal(d$fit[1], 0); expect_equal(d$se[1], 0)
    expect_equal(d$fit, sf$table$A[k] * (1 - exp(-d$tt / sf$table$tau[k])), tolerance = 1e-10)
    expect_true(all(d$lo <= d$fit & d$fit <= d$hi))
  }
})

test_that("the identification flags fire on a ridge and on tau beyond the recording", {
  sf <- list(window = c(0, 25), table = data.frame(
    cell = c("a", "b", "c"), tau = c(10, 40, 8), corr_A_logtau = c(0.95, 0.5, 0.2)))
  fl <- e$dance_traj_sat_flags(sf)
  expect_length(fl, 2L)
  expect_match(fl[1], "In a, A and tau are correlated beyond 0.9")
  expect_match(fl[2], "In b, tau exceeds the 25 h recording")
})

test_that("two-stage: A, tau and A/tau per group are means of the participants' own", {
  set.seed(5)
  p <- data.frame(subject = 1:20, group = factor(rep(c("x", "y"), each = 10)),
                  A_sat = rnorm(20, 20, 2), tau = rgamma(20, 20, 2))
  mod <- list(trend_type = "exp_sat", time_vec = 0:25)
  sm <- e$dance_ts_sat_summary(p, mod)
  expect_true(isTRUE(sm$ok))
  expect_identical(sm$quantities, c("A_sat", "tau", "init_slope"))
  r <- sm$table[sm$table$quantity == "init_slope" & sm$table$cell == "x", ]
  expect_equal(r$estimate, mean(p$A_sat[1:10] / p$tau[1:10]))
  cv <- e$dance_ts_sat_curves(p, mod, n_time = 6)
  d <- cv$table[cv$table$cell == "y", ]
  expect_equal(d$fit, colMeans(t(vapply(11:20, function(i)
    p$A_sat[i] * (1 - exp(-d$tt / p$tau[i])), numeric(6)))), tolerance = 1e-12)
  expect_false(isTRUE(e$dance_ts_sat_summary(p, list(trend_type = "linear"))$ok))
})

test_that("the tau profile uses one random structure and refines its maximum", {
  skip_heavy()
  d <- mk_long(taus = c(A = 9, B = 9), A = c(A = 20, B = 24), nper = 10)
  pr <- e$dance_traj_profile_tau(d, 24, 1, NULL,
                                 tau_grid = exp(seq(log(2), log(60), length.out = 10)))
  expect_true(isTRUE(pr$ok))
  expect_equal(length(unique(pr$profile$rung)), 1L)
  expect_false(pr$flat)
  # the refined maximum is a point of the profile, the best of them, inside the interval
  expect_true(any(abs(pr$profile$tau - pr$tau) < 1e-9))
  expect_equal(max(pr$profile$logLik), pr$profile$logLik[which.min(abs(pr$profile$tau - pr$tau))])
  expect_true(pr$ci[1] < pr$tau && pr$tau < pr$ci[2])
  expect_true(pr$ci[1] < 9 && 9 < pr$ci[2])
  # the fit handed back is at the estimate, and says tau was estimated
  expect_equal(pr$fit$spec$tau, pr$tau)
  expect_true(isTRUE(pr$fit$tau_estimated))
  # a light fit gives the same likelihood as a full one
  sp <- e$dance_traj_spec(d, 24, 1, "exp_sat", tau = 9)
  a <- e$dance_traj_fit(sp, REML = FALSE, light = TRUE)
  b <- e$dance_traj_fit(sp, REML = FALSE)
  expect_equal(as.numeric(stats::logLik(a$model)), as.numeric(stats::logLik(b$model)), tolerance = 1e-6)
})

test_that("the change's factor is spelled out per trend", {
  tw <- e$dance_traj_trend_window("exp_sat", 8, 39, 18.3)
  expect_equal(e$dance_traj_change_factor_text(tw),
               "0.8162 = 1 - e^(-31 h / 18.3 h), the share of A the rise has reached by the last observation")
  expect_match(e$dance_traj_change_factor_text(e$dance_traj_trend_window("linear", 8, 33)),
               "^25, the 25 hours from the first to the last observation$")
  expect_match(e$dance_traj_change_factor_text(e$dance_traj_trend_window("log", 8, 33)),
               "log\\(1 \\+ 25 h\\)")
})

test_that("a search range bounds every per-cell tau and its interval", {
  skip_heavy()
  ff <- mk_fit(mk_long()); skip_if_not(isTRUE(ff$ok))
  free <- e$dance_traj_sat_fit(ff)
  sf <- e$dance_traj_sat_fit(ff, tau_range = c(3, 30))
  skip_if_not(isTRUE(free$ok) && isTRUE(sf$ok))
  expect_identical(sf$tau_map$par, "et")
  st <- sf$table
  expect_true(all(st$tau_lo >= 3 & st$tau_hi <= 30))
  # with the optimum inside the range, the constrained fit IS the free fit
  expect_equal(st$tau, free$table$tau, tolerance = 1e-3)
  expect_equal(st$A, free$table$A, tolerance = 1e-3)
  # A/tau through the logistic map, against a numerical derivative
  V <- e$dance_traj_sat_V(sf$fit)
  nm <- c(e$dance_traj_sat_coef(sf, "A", 2), e$dance_traj_sat_coef(sf, "tau", 2))
  th <- nlme::fixef(sf$fit)[nm]
  g <- function(p) p[1] / sf$tau_map$tau(p[2]); h <- 1e-6
  J <- c((g(th + c(h, 0)) - g(th - c(h, 0))) / (2 * h), (g(th + c(0, h)) - g(th - c(0, h))) / (2 * h))
  expect_equal(st$init_slope_se[2], sqrt(drop(t(J) %*% V[nm, nm] %*% J)), tolerance = 1e-6)
  # the curves' band goes through the same map
  cv <- e$dance_traj_sat_curves(sf, n_time = 5)
  d <- cv$table[cv$table$cell == sf$cells[2], ]
  gc <- function(p, t) p[1] * (1 - exp(-t / sf$tau_map$tau(p[2])))
  tk <- d$tt[3]
  Jc <- c((gc(th + c(h, 0), tk) - gc(th - c(h, 0), tk)) / (2 * h),
          (gc(th + c(0, h), tk) - gc(th - c(0, h), tk)) / (2 * h))
  expect_equal(d$se[3], sqrt(drop(t(Jc) %*% V[nm, nm] %*% Jc)), tolerance = 1e-6)
})

test_that("a range that excludes a cell's tau fails with the reason, not a guess", {
  skip_heavy()
  ff <- mk_fit(mk_long()); skip_if_not(isTRUE(ff$ok))
  sf <- e$dance_traj_sat_fit(ff, tau_range = c(8, 40))   # cell A's planted tau is 5 h
  if (isTRUE(sf$ok)) {
    # if it converged, the estimate sits at the bound and is flagged as one
    expect_true(any(grepl("end of the search range", e$dance_traj_sat_flags(sf), fixed = TRUE)))
  } else {
    expect_match(sf$message, "outside 8 to 40 h", fixed = TRUE)
    expect_match(sf$message, "in A", fixed = TRUE)
  }
})
