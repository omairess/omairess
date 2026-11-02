# ==============================================================================
# tests/testthat/test-model-comparison.R
#
# The model comparison covers every structure -- intercept only, trend only,
# rhythm only, both -- under both approaches, and counts parameters correctly:
# a held tau is not a parameter, a profiled or per-participant tau is. On data
# simulated from a known structure, the comparison should prefer that structure.
# ==============================================================================
.libPaths(c("~/Rlib", .libPaths()))
app_dir <- if (dir.exists("server")) "." else if (dir.exists("../../server")) "../.." else "DANCE"
`%||%` <- function(a, b) if (is.null(a)) b else a

e <- new.env()
for (f in c("07_helpers_circular", "08_helpers_cosinor", "08b_helpers_popcosinor", "08c_helpers_circstat",
            "08d_helpers_traj", "08e_helpers_trajfit", "08f_helpers_trajinf", "08g_helpers_trajband",
            "08h_helpers_twostage", "08i_helpers_trajsat"))
  sys.source(file.path(app_dir, "server", paste0(f, ".R")), envir = e)
ex <- parse(file.path(app_dir, "server", "72_harmonic.R"), encoding = "UTF-8")
for (x in ex) if (is.call(x) && identical(x[[1]], as.name("<-")) && is.name(x[[2]]) &&
                  as.character(x[[2]]) %in% c("fit_cosinor", "fit_cosinor_nonlinear")) eval(x, e)

tp <- seq(0, 30, by = 2)
sim <- function(trend, rhythm, seed) {
  set.seed(seed); g <- rep(c("A", "B"), each = 12)
  Y <- t(sapply(seq_along(g), function(i) 20 + rnorm(1) +
    (if (trend) (3 + 2 * (g[i] == "B")) * (1 - exp(-tp / 10)) else 0) +
    (if (rhythm) 3 * cos(2 * pi * (tp - 8) / 24) else 0) + rnorm(length(tp), 0, 0.5)))
  list(Y = Y, g = g)
}

test_that("the candidate space is every structure, made explicit", {
  ms <- e$dance_cosinor_model_set(2)
  lab <- vapply(ms, function(m) m$label, "")
  st <- vapply(ms, function(m) m$structure, "")
  expect_true(all(c("intercept only", "linear only", "log only", "exp_sat only",
                    "none + H1", "none + H1-H2", "exp_sat + H1", "exp_sat + H1-H2") %in% lab))
  expect_identical(sort(unique(st)), c("combined", "intercept", "rhythm", "trend"))
  expect_length(ms, 12)
})

test_that("two-stage: the comparison prefers the generating structure", {
  for (case in list(list(tr = TRUE, rh = FALSE, best = "trend"),
                    list(tr = FALSE, rh = TRUE, best = "rhythm"),
                    list(tr = TRUE, rh = TRUE, best = "combined"))) {
    s <- sim(case$tr, case$rh, 11)
    tab <- e$dance_ts_model_set(s$Y, tp, 24, 1, "exp_sat", tau_min = 1, tau_max = 60, tau_fixed = 10, k_max = 2)
    expect_identical(tab$structure[1], case$best, info = paste(case$best, tab$model[1]))
    expect_setequal(attr(tab, "candidates"), vapply(e$dance_cosinor_model_set(2), `[[`, "", "label"))
  }
})

test_that("a held tau is not a parameter; an estimated one is", {
  s <- sim(TRUE, TRUE, 12)
  y <- s$Y[1, ]
  f_fix <- e$fit_cosinor(tp, y, 24, 1, "exp_sat", tau_fixed = 10, t_origin = 0)
  f_est <- e$fit_cosinor(tp, y, 24, 1, "exp_sat", tau_min = 1, tau_max = 60, t_origin = 0)
  # intercept + A + 2 harmonic coefficients + residual variance, and tau when estimated
  k_of <- function(f) (f$aic - (-2 * e$dance_gaussian_loglik(f$ss_resid, length(y), stats::sd(y)))) / 2
  expect_equal(k_of(f_fix), 5)
  expect_equal(k_of(f_est), 6)
})

test_that("mixed: ML fits, correct counts, and the generating structure preferred", {
  skip_if_not_installed("lme4")
  for (case in list(list(tr = TRUE, rh = FALSE, best = "trend"),
                    list(tr = TRUE, rh = TRUE, best = "combined"))) {
    s <- sim(case$tr, case$rh, 21)
    d <- e$dance_traj_long(s$Y, tp, sprintf("S%02d", seq_along(s$g)), list(Group = s$g))
    tab <- suppressWarnings(e$dance_traj_model_set(d, 24, k_max = 1, t0 = 0, tau_mode = "fixed",
                                                   tau_value = 10))
    expect_false(is.null(tab))
    expect_identical(tab$structure[1], case$best, info = paste(case$best, tab$model[1]))
    # intercept only: two cell means, a participant variance and the residual
    expect_equal(tab$k[tab$model == "intercept only"], 4)
    # the log-likelihoods are ML: refitting the reported candidate by ML reproduces it
    sp <- e$dance_traj_spec(d, 24, 0, "linear", t0 = 0)
    ft <- suppressWarnings(suppressMessages(e$dance_traj_fit(sp, REML = FALSE)))
    expect_equal(tab$logLik[tab$model == "linear only"], as.numeric(stats::logLik(ft$model)),
                 tolerance = 1e-6)
  }
  # a profiled tau counts as one parameter more than the model's own count
  s <- sim(TRUE, FALSE, 22)
  d <- e$dance_traj_long(s$Y, tp, sprintf("S%02d", seq_along(s$g)), list(Group = s$g))
  tp_est <- suppressWarnings(e$dance_traj_model_set(d, 24, k_max = 1, t0 = 0, tau_mode = "estimate",
                                                    tau_lo = 2, tau_hi = 60))
  r <- tp_est[tp_est$model == "exp_sat only", ]
  expect_true(isTRUE(r$tau_estimated))
  expect_equal(r$k, tp_est$k[tp_est$model == "linear only"] + 1)
})

test_that("mixed: with AR(1) residuals every candidate is fitted at the reported phi", {
  skip_if_not_installed("lme4")
  s <- sim(TRUE, TRUE, 23)
  d <- e$dance_traj_long(s$Y, tp, sprintf("S%02d", seq_along(s$g)), list(Group = s$g))
  t0 <- suppressWarnings(e$dance_traj_model_set(d, 24, k_max = 1, t0 = 0, tau_mode = "fixed", tau_value = 10))
  t1 <- suppressWarnings(e$dance_traj_model_set(d, 24, k_max = 1, t0 = 0, tau_mode = "fixed", tau_value = 10,
                                                phi = 0.5))
  expect_equal(attr(t1, "phi"), 0.5)
  expect_match(attr(t1, "residual_note"), "AR(1)", fixed = TRUE)
  expect_match(attr(t0, "residual_note"), "independent", fixed = TRUE)
  # phi is one more parameter in every candidate
  m <- intersect(t0$model, t1$model)
  expect_equal(t1$k[match(m, t1$model)], t0$k[match(m, t0$model)] + 1)
  # the log-likelihood is the data's: the whitened fit's plus log|det W|
  sp <- e$dance_traj_spec(d, 24, 0, "linear", t0 = 0)
  ft <- suppressWarnings(suppressMessages(e$dance_traj_fit(sp, REML = FALSE, phi = 0.5)))
  expect_equal(t1$logLik[t1$model == "linear only"],
               as.numeric(stats::logLik(ft$model)) + e$dance_car1_logdet(sp, 0.5), tolerance = 1e-6)
})
