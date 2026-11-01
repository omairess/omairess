# ==============================================================================
# tests/testthat/test-outcomes-calculus.R
#
# The three-domain outcomes are defined once (server/08_helpers_cosinor.R) and
# read by both estimators. These pin the definitions to independent
# computations -- numerical integration and differentiation of the trend and of
# the fitted curve -- and the relations the brief states: a straight line's
# initial and average rates are one number; with a common tau and window, A,
# the initial rate and the total change are proportional.
# ==============================================================================
.libPaths(c("~/Rlib", .libPaths()))
app_dir <- if (dir.exists("server")) "." else if (dir.exists("../../server")) "../.." else "DANCE"
`%||%` <- function(a, b) if (is.null(a)) b else a

e <- new.env()
for (f in c("07_helpers_circular", "08_helpers_cosinor", "08c_helpers_circstat", "08d_helpers_traj", "08e_helpers_trajfit",
            "08f_helpers_trajinf", "08g_helpers_trajband", "08h_helpers_twostage"))
  sys.source(file.path(app_dir, "server", paste0(f, ".R")), envir = e)
ex <- parse(file.path(app_dir, "server", "72_harmonic.R"), encoding = "UTF-8")
for (x in ex) if (is.call(x) && identical(x[[1]], as.name("<-")) && is.name(x[[2]]) &&
                  as.character(x[[2]]) %in% c("fit_cosinor", "fit_cosinor_nonlinear")) eval(x, e)

test_that("the trend factors are the integral, the endpoint and the slope of the shape", {
  for (tr in c("linear", "log", "exp_sat")) for (T in c(4, 22, 44)) {
    tau <- if (tr == "exp_sat") 13 else NULL
    f <- e$dance_trend_factors(tr, T, tau)
    s <- function(t) e$dance_trend_shape(tr, t, tau)
    expect_equal(f$end, s(T) - s(0), tolerance = 1e-12, info = tr)
    expect_equal(f$mean, stats::integrate(s, 0, T, rel.tol = 1e-12)$value / T,
                 tolerance = 1e-9, info = paste(tr, T))
    h <- 1e-6
    expect_equal(f$rate0, (s(h) - s(0)) / h, tolerance = 1e-5, info = tr)
    expect_equal(f$avg_rate, f$end / T, tolerance = 1e-12, info = tr)
  }
  expect_null(e$dance_trend_factors("none", 10))
  expect_null(e$dance_trend_factors("exp_sat", 10, tau = NA))
})

test_that("a straight line's initial and average rates are one number", {
  o <- e$dance_trend_outcomes("linear", 0.37, T = 31)
  expect_equal(o$rate0, o$avg_rate)
  expect_equal(o$change, 0.37 * 31)
  expect_true(is.na(o$asymptote) && is.na(o$t50))
})

test_that("with a common tau and window, A, the initial rate and the change are proportional", {
  T <- 31; tau <- 18
  A <- c(2, 5, 9.5)
  o <- lapply(A, function(a) e$dance_trend_outcomes("exp_sat", a, T, tau))
  ch <- vapply(o, `[[`, 0, "change"); r0 <- vapply(o, `[[`, 0, "rate0")
  expect_equal(ch / A, rep(1 - exp(-T / tau), 3))
  expect_equal(r0 / A, rep(1 / tau, 3))
  expect_equal(o[[1]]$t50, tau * log(2))
  expect_equal(o[[1]]$frac_reached, 1 - exp(-T / tau))
})

test_that("the window mean of a rhythm is its integral, and zero over whole cycles", {
  bc <- c(2, -0.5); bs <- c(1, 0.7)
  for (T in c(24, 31, 44)) {
    hm <- e$dance_harmonic_mean_factors(2, 24, T)
    yh <- function(t) bc[1] * cos(2 * pi * t / 24) + bs[1] * sin(2 * pi * t / 24) +
                      bc[2] * cos(4 * pi * t / 24) + bs[2] * sin(4 * pi * t / 24)
    expect_equal(sum(bc * hm$cos + bs * hm$sin),
                 stats::integrate(yh, 0, T, rel.tol = 1e-12)$value / T, tolerance = 1e-9)
  }
  hm24 <- e$dance_harmonic_mean_factors(2, 24, 48)
  expect_equal(c(hm24$cos, hm24$sin), rep(0, 4), tolerance = 1e-12)
  expect_length(e$dance_harmonic_mean_factors(0, 24, 30)$cos, 0)
})

test_that("the time-averaged fitted response is the integral of the fitted curve", {
  lo <- e$dance_level_outcomes(20, "exp_sat", 6, tau = 12, beta_cos = 3, beta_sin = -1,
                               period = 24, T = 31)
  yhat <- function(t) 20 + 6 * (1 - exp(-t / 12)) + 3 * cos(2 * pi * t / 24) - sin(2 * pi * t / 24)
  expect_equal(lo$fitted_mean, stats::integrate(yhat, 0, 31, rel.tol = 1e-12)$value / 31,
               tolerance = 1e-9)
  expect_equal(lo$initial_midline, 20)
  expect_equal(lo$mean_midline, 20 + 6 * (1 - (12 / 31) * (1 - exp(-31 / 12))))
})

test_that("the model structure is read from the trend and the harmonic count", {
  expect_identical(e$dance_model_structure("none", 2), "rhythm")
  expect_identical(e$dance_model_structure("log", 0), "trend")
  expect_identical(e$dance_model_structure("exp_sat", 1), "combined")
  expect_identical(e$dance_model_structure("none", 0), "intercept")
})

# ---- the mixed model: every structure, every trend -------------------------
set.seed(2)
tp <- seq(8, 38, by = 2); grp <- rep(c("A", "B"), each = 12)
Y <- t(sapply(seq_along(grp), function(i) 20 + rnorm(1) + (3 + 2 * (grp[i] == "B")) *
                (1 - exp(-(tp - 8) / 10)) + 3 * cos(2 * pi * (tp - 16) / 24) + rnorm(length(tp), 0, 0.5)))
d <- e$dance_traj_long(Y, tp, sprintf("S%02d", seq_along(grp)), list(Group = grp))

test_that("mixed: the cell outcomes are the closed forms of the fitted cell curve", {
  skip_if_not_installed("lme4")
  for (tr in c("none", "linear", "log", "exp_sat")) for (K in 0:2) {
    if (tr == "none" && K == 0) next
    sp <- e$dance_traj_spec(d, 24, K, tr, tau = 10)
    ff <- suppressWarnings(suppressMessages(e$dance_traj_fit(sp)))
    expect_true(isTRUE(ff$ok), info = paste(tr, K))
    ls <- e$dance_traj_level_summary(ff)
    q <- unique(ls$table$quantity)
    # each quantity appears only where it is distinct
    expect_identical("fitted_mean" %in% q, K > 0, info = paste(tr, K))
    expect_identical(all(c("mesor", "change", "avg_rate") %in% q), tr != "none", info = paste(tr, K))
    expect_identical("rate0" %in% q, tr %in% c("log", "exp_sat"), info = paste(tr, K))
    bv <- e$dance_traj_beta(ff); grid <- e$dance_traj_cell_grid(ff$spec)
    tg <- seq(8, 38, length.out = 20001)
    for (i in seq_len(nrow(grid))) {
      X <- e$dance_traj_design_rows(ff, grid[i, , drop = FALSE], tg)
      yh <- drop(X[, names(bv$beta)] %*% bv$beta)
      est <- function(nm) ls$table$estimate[ls$table$quantity == nm & ls$table$cell == grid$.cell[i]]
      expect_equal(if (K > 0) est("fitted_mean") else if (tr != "none") est("mesor") else NA,
                   mean((head(yh, -1) + tail(yh, -1)) / 2), tolerance = 1e-6, info = paste(tr, K))
      if (tr != "none") {
        expect_equal(est("avg_rate") * 30, est("change"), tolerance = 1e-10, info = paste(tr, K))
        if (tr != "linear") {
          r0 <- est("rate0"); ch <- est("change")
          f <- e$dance_trend_factors(tr, 30, if (tr == "exp_sat") 10)
          expect_equal(r0 / f$rate0, ch / f$end, tolerance = 1e-8, info = paste(tr, K))
        }
      }
    }
    # a trend-only model offers no rhythm comparison
    ws <- e$dance_traj_pair_whats(ff)
    if (K == 0) expect_false(any(grepl("^(amplitude|phase)|circadian|shape", ws)), info = tr)
  }
})

test_that("mixed: the change, average-rate and initial-rate tests are one test", {
  skip_if_not_installed("lme4")
  sp <- e$dance_traj_spec(d, 24, 1, "exp_sat", tau = 10)
  ff <- suppressWarnings(suppressMessages(e$dance_traj_fit(sp)))
  ct <- lapply(c("change", "avg_rate", "rate0", "trend"), function(w)
    e$dance_traj_contrasts(ff, w, adjust = "none", df_method = "satterthwaite"))
  for (k in 2:4) expect_equal(ct[[k]]$table$p_raw, ct[[1]]$table$p_raw, tolerance = 1e-8)
})

# ---- two-stage: every structure, every trend --------------------------------
test_that("two-stage: per-participant outcomes for every structure, missing values included", {
  set.seed(3); tv <- seq(0, 30, by = 2); n <- 20; g <- rep(c("A", "B"), each = 10)
  Y2 <- t(sapply(1:n, function(i) 20 + rnorm(1) + (3 + 2 * (g[i] == "B")) * (1 - exp(-tv / 12)) +
                   2 * cos(2 * pi * (tv - 6) / 24) + rnorm(length(tv), 0, .5)))
  Y2[3, 5] <- NA; Y2[7, c(1, 9)] <- NA
  for (tr in c("none", "linear", "log", "exp_sat")) for (K in 0:1) {
    if (tr == "none" && K == 0) next
    fa <- e$dance_ts_fit_all(Y2, tv, 24, K, tr, tau_min = 1, tau_max = 60)
    ok <- vapply(fa$fits, function(f) isTRUE(f$success), TRUE)
    expect_true(all(ok), info = paste(tr, K))
    mod <- list(individual_params = fa$params, individual_fits = fa$fits, period = 24,
                n_harmonics = K, trend_type = tr, time_vec = tv, t_offset = 0,
                tau_spec = list(mode = "estimate"))
    df <- e$dance_ts_frame(mod, g)
    ls <- e$dance_ts_level_summary(df, mod)
    expect_identical(unique(ls$table$quantity), names(e$dance_ts_level_columns(mod)), info = paste(tr, K))
    if (tr != "none") expect_equal(df$avg_rate * 30, df$trend_change, tolerance = 1e-8, info = paste(tr, K))
    # each participant's fitted mean is the integral of their own fitted curve
    i <- 7; fi <- fa$fits[[i]]
    yh <- function(t) e$dance_rhythm_from_coefs(fi$coefs, t, 24, K, tr, include_trend = TRUE,
                                                t_offset = fi$t_offset %||% 0)
    tg <- seq(0, 30, length.out = 6001)
    yv <- vapply(tg, yh, 0)
    expect_equal(df$fitted_mean[i], mean((head(yv, -1) + tail(yv, -1)) / 2), tolerance = 1e-5,
                 info = paste(tr, K))
    for (w in e$dance_ts_pair_whats(mod)) {
      pr <- e$dance_ts_pairwise(df, mod, w)
      expect_true(isTRUE(pr$ok), info = paste(tr, K, w))
    }
  }
})
