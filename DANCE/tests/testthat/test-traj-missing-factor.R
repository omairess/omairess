# ==============================================================================
# tests/testthat/test-traj-missing-factor.R
#
# A participant with no value for a design factor cannot be placed in a cell.
# lme4 dropped those rows on its own, while the frame everything else is built
# from still carried them: the residuals were longer than the fitted values,
# and a participant curve fell in a cell "NA" that has no colour, which crashed
# the Fitted curves plot. dance_traj_spec() now drops them once and says so.
# Also: the R-squared keeps Z sparse (a dense Z was gigabytes on large data).
# ==============================================================================
.libPaths(c("~/Rlib", .libPaths()))
app_dir <- if (dir.exists("server")) "." else if (dir.exists("../../server")) "../.." else "DANCE"
`%||%` <- function(a, b) if (is.null(a)) b else a
e <- new.env()
for (f in c("08_helpers_cosinor", "08c_helpers_circstat", "08d_helpers_traj",
            "08e_helpers_trajfit", "08f_helpers_trajinf", "08g_helpers_trajband"))
  sys.source(file.path(app_dir, "server", paste0(f, ".R")), envir = e)

mk <- function(n = 24, seed = 4) {
  set.seed(seed)
  tp <- seq(0, 22, by = 2)
  g <- rep(c("A", "B"), length.out = n); g[c(3, 8)] <- NA
  Y <- t(sapply(seq_len(n), function(i) 30 + rnorm(1, 0, 2) + 3 * (isTRUE(g[i] == "B")) +
                  6 * cos(2 * pi * (tp - 4) / 24) + rnorm(length(tp), 0, 1)))
  e$dance_traj_long(Y, tp, sprintf("S%02d", seq_len(n)), list(Group = g))
}

test_that("participants with no group are excluded once, and the exclusion is reported", {
  skip_if_not_installed("lme4")
  d <- mk()
  sp <- e$dance_traj_spec(d, 24, 1, "none")
  expect_true(isTRUE(sp$ok))
  expect_false(any(is.na(sp$data$Group)))
  expect_equal(nlevels(droplevels(sp$data$subject)), 22L)
  expect_match(sp$missing_note, "2 participants with no value for 'Group' were excluded")
  expect_match(sp$missing_note, "S03, S08")
  ff <- suppressWarnings(suppressMessages(e$dance_traj_fit(sp)))
  expect_true(isTRUE(ff$ok))
  # the frame and the model agree, row for row
  r <- e$dance_traj_residuals(ff)
  expect_equal(length(r$raw), nrow(sp$data))
  expect_equal(length(r$fitted), stats::nobs(ff$model))
  # every participant curve lies in a cell the model has
  pp <- e$dance_traj_participant_predict(ff, n_time = 10)
  expect_true(all(pp$table$cell %in% c("A", "B")))
  # the R-squared is unchanged by keeping Z sparse: compare with the dense sum
  r2 <- e$dance_traj_r2(ff)
  Z <- as.matrix(lme4::getME(ff$model, "Z")); vc <- lme4::VarCorr(ff$model)
  S <- as.matrix(vc[[1]]); q <- ncol(S); nl <- ncol(Z) / q
  vr <- 0; for (l in seq_len(nl)) { Zl <- Z[, (l - 1) * q + seq_len(q), drop = FALSE]; vr <- vr + rowSums((Zl %*% S) * Zl) }
  expect_equal(r2$var_random, mean(vr), tolerance = 1e-10)
})

test_that("a frame with no missing group carries no note", {
  skip_if_not_installed("lme4")
  d <- mk(); d <- d[!is.na(d$Group), ]
  sp <- e$dance_traj_spec(d, 24, 1, "none")
  expect_null(sp$missing_note)
})

test_that("a blank group label is a missing value, not a group called \"\"", {
  skip_if_not_installed("lme4")
  set.seed(2); tp <- seq(0, 22, by = 2)
  g <- rep(c("A", "B"), 10); g[c(1, 6)] <- c("", "  ")
  Y <- t(sapply(seq_along(g), function(i) 30 + rnorm(1) + 6 * cos(2 * pi * (tp - 4) / 24) + rnorm(length(tp))))
  d <- e$dance_traj_long(Y, tp, sprintf("S%02d", seq_along(g)), list(Group = g))
  expect_true(all(is.na(d$Group[d$subject %in% c("S01", "S06")])))
  expect_setequal(levels(d$Group), c("A", "B"))
  sp <- e$dance_traj_spec(d, 24, 1, "none")
  expect_match(sp$missing_note, "2 participants")
})

test_that("two-stage: a blank group label puts the participant in no group", {
  e2 <- new.env()
  sys.source(file.path(app_dir, "server", "08h_helpers_twostage.R"), envir = e2)
  mod <- list(individual_params = data.frame(subject = 1:4, mesor = 1:4, trend_change = 0),
              trend_type = "none")
  df <- e2$dance_ts_frame(mod, c("A", "", "B", "A"))
  expect_equal(nrow(df), 3L)
  expect_setequal(levels(df$group), c("A", "B"))
})
