# ==============================================================================
# tests/testthat/test-twostage-helpers.R
#
# The two-stage comparison kernels in server/08h_helpers_twostage.R, on planted
# data: the per-curve frame, the MANOVA primary, the component tests, the
# pairwise tests and the group curves with their SE band. Nothing here is a new
# procedure -- stats::manova, one-way ANOVA, Bingham's population-mean cosinor,
# Watson-Williams, Welch's t -- so the checks are that the right kernel is
# called on the right column and that a planted difference is found while a
# planted equality is not.
# ==============================================================================
.libPaths(c("~/Rlib", .libPaths()))
app_dir <- if (dir.exists("server")) "." else if (dir.exists("../../server")) "../.." else "DANCE"
`%||%` <- function(a, b) if (is.null(a)) b else a
e <- new.env()
for (f in c("07_helpers_circular", "08_helpers_cosinor", "08b_helpers_popcosinor",
            "08c_helpers_circstat", "08d_helpers_traj", "08e_helpers_trajfit",
            "08f_helpers_trajinf", "08g_helpers_trajband", "08h_helpers_twostage"))
  sys.source(file.path(app_dir, "server", paste0(f, ".R")), envir = e)

# a two-stage model object as the run button stores it: per-participant OLS
# parameters, the per-participant fits with their coefficient vectors, and the
# group vector means
mk_mod <- function(n = c(A = 24, B = 18), acro = c(A = 15, B = 18), amp = c(A = 8, B = 8),
                   mesor = c(A = 50, B = 50), seed = 21, trend = "none") {
  set.seed(seed)
  g <- rep(names(n), n); N <- sum(n)
  phi <- 2 * pi * (acro[g] + rnorm(N, 0, 0.5)) / 24
  A <- amp[g] + rnorm(N, 0, 0.7); m0 <- mesor[g] + rnorm(N, 0, 1.5)
  slope <- if (trend == "linear") rnorm(N, 0.1, 0.02) else numeric(0)
  a <- A * cos(phi); b <- A * sin(phi)
  ip <- data.frame(subject = seq_len(N), mesor = m0, beta_cos_1 = a, beta_sin_1 = b,
                   amplitude_1 = A, acrophase_rad_1 = phi %% (2 * pi),
                   acrophase_time_1 = (phi %% (2 * pi)) * 24 / (2 * pi),
                   r_squared = runif(N, .7, .95), p_value = 1e-4, stringsAsFactors = FALSE)
  if (trend == "linear") ip$trend_linear <- slope
  fits <- lapply(seq_len(N), function(i) list(
    success = TRUE, t_offset = 0,
    coefs = if (trend == "linear") c(m0[i], slope[i], a[i], b[i]) else c(m0[i], a[i], b[i])))
  gf <- lapply(names(n), function(lv) {
    i <- g == lv
    list(n = sum(i), mean_coefs = if (trend == "linear")
      c(mean(m0[i]), mean(slope[i]), mean(a[i]), mean(b[i]))
      else c(mean(m0[i]), mean(a[i]), mean(b[i])))
  })
  names(gf) <- names(n)
  list(mod = list(approach = "two_stage", period = 24, n_harmonics = 1L, trend_type = trend,
                  individual_params = ip, individual_fits = fits, group_fits = gf,
                  time_vec = 0:23, origin_shift = 0, t_offset = 0),
       g = factor(g))
}

test_that("the frame attaches the group by participant and drops unlabelled rows", {
  m <- mk_mod()
  df <- e$dance_ts_frame(m$mod, m$g)
  expect_equal(nrow(df), 42L)
  expect_equal(as.character(df$group), as.character(m$g))
  gg <- as.character(m$g); gg[1:3] <- NA
  df2 <- e$dance_ts_frame(m$mod, gg)
  expect_equal(nrow(df2), 39L)
  # no grouping at all: one group, nothing dropped
  df3 <- e$dance_ts_frame(m$mod, NULL)
  expect_equal(nlevels(df3$group), 1L)
  expect_equal(e$dance_ts_coef_cols(m$mod), c("mesor", "beta_cos_1", "beta_sin_1"))
  expect_equal(e$dance_ts_coef_cols(mk_mod(trend = "linear")$mod),
               c("mesor", "trend_linear", "beta_cos_1", "beta_sin_1"))
})

test_that("the MANOVA primary is Wilks' lambda from stats::manova and finds a planted phase shift", {
  m <- mk_mod(); df <- e$dance_ts_frame(m$mod, m$g)
  r <- e$dance_ts_manova(df, e$dance_ts_coef_cols(m$mod))
  expect_true(isTRUE(r$ok))
  expect_match(r$method, "Wilks")
  expect_equal(r$n_par, 3L); expect_equal(r$n_groups, 2L)
  expect_lt(r$p, 1e-4)
  # the same numbers stats::manova gives directly
  Y <- as.matrix(df[, c("mesor", "beta_cos_1", "beta_sin_1")])
  st <- summary(stats::manova(Y ~ df$group), test = "Wilks")$stats
  expect_equal(r$F, unname(st[1, "approx F"]))
  expect_equal(r$wilks, unname(st[1, "Wilks"]))
  # groups smaller than the coefficient count are refused, not fitted
  small <- df[c(1:3, 25:27), ]
  expect_false(isTRUE(e$dance_ts_manova(small, e$dance_ts_coef_cols(m$mod))$ok))
})

test_that("the component tests come from the named kernels and answer the planted design", {
  m <- mk_mod(acro = c(A = 15, B = 18), amp = c(A = 8, B = 8), mesor = c(A = 50, B = 50))
  df <- e$dance_ts_frame(m$mod, m$g)
  rows <- e$dance_ts_components(df, m$mod)
  lab <- vapply(rows, function(r) r$label, character(1))
  expect_true("Constant term (beta_0)" %in% lab)
  expect_true("H1 amplitude" %in% lab); expect_true("H1 acrophase" %in% lab)
  expect_true(any(grepl("Watson-Williams", lab)))
  get <- function(l) rows[[match(l, lab)]]
  expect_equal(get("Constant term (beta_0)")$method, "one-way ANOVA")
  expect_match(get("H1 amplitude")$method, "Bingham")
  expect_equal(get("H1 acrophase (Watson-Williams, unweighted)")$method, "Watson-Williams F")
  # equal constants: not rejected; different phases: rejected on both phase tests
  expect_gt(get("Constant term (beta_0)")$p, 0.05)
  expect_lt(get("H1 acrophase")$p, 1e-3)
  expect_lt(get("H1 acrophase (Watson-Williams, unweighted)")$p, 1e-3)
  # Bingham's caution: with the phases differing, the amplitude row is flagged
  expect_false(is.null(get("H1 amplitude")$note))
  # and it is NOT flagged when the phases agree and the amplitude differs
  m2 <- mk_mod(acro = c(A = 15, B = 15), amp = c(A = 8, B = 11))
  rows2 <- e$dance_ts_components(e$dance_ts_frame(m2$mod, m2$g), m2$mod)
  lab2 <- vapply(rows2, function(r) r$label, character(1))
  amp2 <- rows2[[match("H1 amplitude", lab2)]]
  expect_null(amp2$note); expect_lt(amp2$p, 1e-3)
  # the ANOVA on the constant is the one-way F from the app's own kernel
  lt <- e$dance_group_linear_test(df$mesor, df$group)
  expect_equal(get("Constant term (beta_0)")$F, lt$F)
})

test_that("pairwise: Welch's t for scalars, Watson-Williams for phase, Holm across the family", {
  m <- mk_mod(n = c(A = 20, B = 16, C = 14), acro = c(A = 15, B = 18, C = 15),
              amp = c(A = 8, B = 8, C = 11))
  df <- e$dance_ts_frame(m$mod, m$g)
  whats <- e$dance_ts_pair_whats(m$mod)
  expect_true(all(c("level", "amplitude1", "phase1", "r_squared") %in% whats))
  pa <- e$dance_ts_pairwise(df, m$mod, "amplitude1", adjust = "holm")
  expect_true(isTRUE(pa$ok)); expect_false(pa$circular)
  expect_equal(nrow(pa$table), 3L)
  expect_equal(pa$table$p_adj, stats::p.adjust(pa$table$p_raw, "holm"))
  # A vs B equal amplitude, C larger than both
  ab <- pa$table[pa$table$cell1 == "A" & pa$table$cell2 == "B", ]
  ac <- pa$table[pa$table$cell1 == "A" & pa$table$cell2 == "C", ]
  expect_gt(ab$p_adj, 0.05); expect_lt(ac$p_adj, 0.01)
  # the row IS a Welch test
  tt <- stats::t.test(df$amplitude_1[df$group == "A"], df$amplitude_1[df$group == "C"], var.equal = FALSE)
  expect_equal(ac$statistic, unname(tt$statistic)); expect_equal(ac$p_raw, tt$p.value)
  # phase: circular, no interval, hours on the shortest arc
  pp <- e$dance_ts_pairwise(df, m$mod, "phase1", adjust = "none")
  expect_true(pp$circular); expect_true(all(is.na(pp$table$lo)))
  ab_p <- pp$table[pp$table$cell1 == "A" & pp$table$cell2 == "B", ]
  expect_lt(abs(ab_p$estimate - (-3)), 0.6)       # 15:00 vs 18:00 -> about -3 h
  expect_lt(ab_p$p_raw, 1e-3)
  ac_p <- pp$table[pp$table$cell1 == "A" & pp$table$cell2 == "C", ]
  expect_gt(ac_p$p_raw, 0.05)
  expect_equal(pp$table$p_adj, pp$table$p_raw)
  # an unknown parameter is refused with a message, not an error
  expect_false(isTRUE(e$dance_ts_pairwise(df, m$mod, "nonsense")$ok))
})

test_that("group curves: the line is the mean-coefficient curve, the band is +/- z * SE(t)", {
  m <- mk_mod(trend = "linear")
  times <- seq(0, 23, by = 0.5)
  gc <- e$dance_ts_group_curves(m$mod, m$g, times, include_trend = TRUE, conf = 0.95)
  expect_true(isTRUE(gc$ok)); expect_setequal(gc$cells, c("A", "B"))
  expect_equal(gc$band, "pointwise")
  for (lv in c("A", "B")) {
    d <- gc$table[gc$table$cell == lv, ]
    expect_equal(nrow(d), length(times))
    # the line: the group's mean coefficients through the same predictor
    line <- e$dance_rhythm_from_coefs(m$mod$group_fits[[lv]]$mean_coefs, times, 24, 1L,
                                      "linear", include_trend = TRUE, t_offset = 0)
    expect_equal(d$fit, line, tolerance = 1e-10)
    # the band: SD across participants' own curves / sqrt(n), times z
    i <- which(m$g == lv)
    M <- t(sapply(i, function(k) e$dance_rhythm_from_coefs(m$mod$individual_fits[[k]]$coefs, times, 24, 1L,
                                                            "linear", include_trend = TRUE, t_offset = 0)))
    se <- apply(M, 2, sd) / sqrt(length(i))
    expect_equal(unname(d$se), unname(se), tolerance = 1e-10)
    expect_equal(unname(d$hi - d$fit), unname(qnorm(0.975) * se), tolerance = 1e-10)
    expect_equal(d$n, rep(length(i), length(times)))
    # with linear coefficients the mean of the curves equals the curve of the means
    expect_equal(d$mean_of_curves, d$fit, tolerance = 1e-10)
  }
  # no grouping: one "(all)" cell whose line is the mean of everyone's curve
  g0 <- e$dance_ts_group_curves(m$mod, NULL, times)
  expect_equal(g0$cells, "(all)")
})
