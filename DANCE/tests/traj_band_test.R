# ==============================================================================
# tests/traj_band_test.R — P21 phase 3: bands, difference curves, contrasts
# ==============================================================================
# Layer E is entirely linear in the fixed effects, so almost everything here can
# be checked against an INDEPENDENT computation rather than against itself: a
# band's SE against a hand-built quadratic form, a difference curve against the
# contrast's own variance, the adapter's coefficients against emmeans.
#
# The checks that matter most are the ones about what is NOT true:
#   - a difference band is not two bands subtracted;
#   - dropping the cross-cell covariance is not conservative;
#   - a simultaneous band is wider than a pointwise one, always;
#   - the adapter does not rename a standard error as a standard deviation.
# ==============================================================================

`%||%` <- function(a, b) if (is.null(a)) b else a
app <- if (dir.exists("server")) "." else ".."
e <- new.env()
for (f in c("08_helpers_cosinor", "08c_helpers_circstat", "08d_helpers_traj",
            "08e_helpers_trajfit", "08f_helpers_trajinf", "08g_helpers_trajband"))
  sys.source(file.path(app, "server", paste0(f, ".R")), envir = e)

pass <- 0L; fail <- 0L
chk <- function(ok, msg_ok, msg_bad) {
  if (isTRUE(ok)) { pass <<- pass + 1L; cat("ok   ", msg_ok, "\n") }
  else { fail <<- fail + 1L; cat("FAIL ", msg_bad, "\n") }
}

make <- function(seed, n_per_group = 8, conds = c("p", "q"), nt = 12, amp = 4,
                 shift = 0, shift_cell = NULL, level = 0, level_group = NULL,
                 noise = 1, slope = 0) {
  set.seed(seed); tp <- seq(0, 22, length.out = nt)
  subj <- sprintf("S%03d", seq_len(n_per_group * 2))
  sg <- rep(c("a", "b"), each = n_per_group)
  rows <- list(); k <- 1L
  meta <- list(subject = character(0), Group = character(0), Condition = character(0))
  for (i in seq_along(subj)) {
    sb <- rnorm(1, 0, 2); sa <- rnorm(1, 0, .4); sp <- rnorm(1, 0, .2)
    for (cc in conds) {
      base <- 10 + sb + rnorm(1, 0, 2); a <- amp + sa + rnorm(1, 0, .4)
      ph <- 2 + sp + rnorm(1, 0, .2)
      if (!is.null(shift_cell) && paste(sg[i], cc) == shift_cell) ph <- ph + shift
      lv <- if (!is.null(level_group) && sg[i] == level_group) level else 0
      rows[[k]] <- base + lv + a * cos(2 * pi * tp / 24 - ph) + slope * tp +
                   rnorm(nt, 0, noise)
      meta$subject <- c(meta$subject, subj[i]); meta$Group <- c(meta$Group, sg[i])
      meta$Condition <- c(meta$Condition, cc); k <- k + 1L
    }
  }
  list(Y = do.call(rbind, rows), t = tp, meta = meta)
}
build <- function(g) e$dance_traj_long(g$Y, g$t, g$meta$subject,
                                       list(Group = g$meta$Group, Condition = g$meta$Condition))

suppressWarnings(suppressMessages({
fit <- e$dance_traj_fit(e$dance_traj_spec(build(make(11)), 24, 1, "none"))
stopifnot(isTRUE(fit$ok))

# ---------------------------------------------------------------- the cell grid
cat("-- the cell grid is the full factorial, in a stable order ----------\n")
g <- e$dance_traj_cell_grid(fit$spec)
chk(nrow(g) == 4 && identical(sort(g$.cell), sort(c("a x p", "a x q", "b x p", "b x q"))),
    "a 2 x 2 design gives 4 cells named by their levels",
    sprintf("the grid is wrong: %s", paste(g$.cell, collapse = ", ")))

# ---------------------------------------------------------------- design rows
cat("\n-- the design row is rebuilt on the fit's own terms -----------------\n")
X <- e$dance_traj_design_rows(fit, g[1, , drop = FALSE], c(0, 6, 12))
chk(nrow(X) == 3 && all(colnames(X) %in% names(e$dance_traj_beta(fit)$beta)) &&
      length(colnames(X)) == length(e$dance_traj_beta(fit)$beta),
    "every column of the rebuilt design row is a column the model has, and none is missing",
    "the rebuilt design row does not match the model's fixed effects")
# at t = t0 the harmonic columns are exactly (1, 0) -- the property the adapter's
# intercept recovery depends on
X0 <- e$dance_traj_design_rows(fit, g[1, , drop = FALSE], fit$spec$t0)
chk(abs(X0[1, "c1"] - 1) < 1e-12 && abs(X0[1, "s1"]) < 1e-12,
    "at t = t0 the harmonic columns are exactly (1, 0)",
    "the basis is not anchored where the spec says it is")

# ---------------------------------------------------------------- prediction
cat("\n-- one trajectory per cell, and the band says which kind it is -----\n")
pw <- e$dance_traj_predict(fit, band = "pointwise", n_time = 50)
sm <- e$dance_traj_predict(fit, band = "simultaneous", n_time = 50)
chk(isTRUE(pw$ok) && length(unique(pw$table$cell)) == 4 && nrow(pw$table) == 200,
    "four cells x 50 times, each with a fit and a band",
    "the prediction table is the wrong shape")
chk(sm$multiplier > pw$multiplier &&
      all(sm$table$hi - sm$table$lo >= pw$table$hi - pw$table$lo - 1e-9),
    sprintf("the simultaneous band is wider everywhere (multiplier %.2f vs %.2f)",
            sm$multiplier, pw$multiplier),
    "a simultaneous band came back no wider than a pointwise one")
chk(grepl("POINTWISE", pw$note) && grepl("not a", pw$note) &&
      grepl("SIMULTANEOUS", sm$note),
    "each band states which kind it is and what it does not license",
    "a band came back without saying which kind it is")

# the SE is checked against an independently built quadratic form
bv <- e$dance_traj_beta(fit)
Xc <- e$dance_traj_design_rows(fit, g[2, , drop = FALSE], pw$times)
kp <- intersect(colnames(Xc), names(bv$beta))
se_manual <- sqrt(diag(Xc[, kp, drop = FALSE] %*% bv$V[kp, kp] %*% t(Xc[, kp, drop = FALSE])))
se_shipped <- pw$table$se[pw$table$cell == g$.cell[2]]
chk(max(abs(se_manual - se_shipped)) < 1e-10,
    sprintf("the band's SE matches an independently built x' V x (max diff %.2e)",
            max(abs(se_manual - se_shipped))),
    "the band's SE does not match the quadratic form it claims to be")

# ---------------------------------------------------------------- difference
cat("\n-- a difference curve is ONE contrast, not two curves subtracted ----\n")
dc <- e$dance_traj_diff_curve(fit, "a x p", "a x q", n_time = 50, band = "pointwise")
chk(isTRUE(dc$ok) && nrow(dc$table) == 50,
    "the difference curve is returned over the requested grid",
    "the difference curve failed")
p1 <- pw$table[pw$table$cell == "a x p", ]; p2 <- pw$table[pw$table$cell == "a x q", ]
chk(max(abs(dc$table$diff - (p1$fit - p2$fit))) < 1e-10,
    "its point estimate IS the difference of the two fitted curves",
    "the difference curve does not equal the difference of the curves")
naive <- sqrt(p1$se^2 + p2$se^2)
chk(max(abs(dc$table$se - naive)) > 1e-6,
    sprintf("but its SE is NOT sqrt(se1^2 + se2^2) -- they differ by up to %.3f",
            max(abs(dc$table$se - naive))),
    "the difference SE equals the independent-errors formula, so the covariance was dropped")
# and the exact contrast variance, built independently
D <- e$dance_traj_design_rows(fit, g[match("a x p", g$.cell), , drop = FALSE], dc$table$t) -
     e$dance_traj_design_rows(fit, g[match("a x q", g$.cell), , drop = FALSE], dc$table$t)
se_d <- sqrt(diag(D[, kp, drop = FALSE] %*% bv$V[kp, kp] %*% t(D[, kp, drop = FALSE])))
chk(max(abs(se_d - dc$table$se)) < 1e-10,
    sprintf("it matches the contrast's own quadratic form (max diff %.2e)",
            max(abs(se_d - dc$table$se))),
    "the difference SE is not the contrast variance it claims to be")
chk(isFALSE(e$dance_traj_diff_curve(fit, "a x p", "a x p")$ok) &&
      isFALSE(e$dance_traj_diff_curve(fit, "a x p", "nonesuch")$ok),
    "a cell against itself, and an unknown cell, are both refused",
    "an impossible difference was accepted")

# a planted shift must show up as a difference band that excludes zero
fit_s <- e$dance_traj_fit(e$dance_traj_spec(
  build(make(12, shift = 1.2, shift_cell = "a q")), 24, 1, "none"))
dcs <- e$dance_traj_diff_curve(fit_s, "a x p", "a x q", n_time = 60)
chk(isTRUE(dcs$any_separation),
    sprintf("a planted 1.2 rad shift separates the simultaneous difference band from zero at %d of %d times",
            sum(dcs$table$excludes_zero), nrow(dcs$table)),
    "a planted shift did not separate the difference band from zero")
dcn <- e$dance_traj_diff_curve(fit, "a x p", "b x p", n_time = 60)
chk(!isTRUE(dcn$any_separation),
    "and with nothing planted the band covers zero throughout",
    "a null difference band excluded zero somewhere")

# ------------------------------------------------- cross-cell covariance
cat("\n-- the cross-cell covariance is exact, and dropping it is not safe -\n")
co <- e$dance_traj_cell_coefs(fit, 1)
pcv <- e$dance_traj_pair_cov(co, 1, 2)
chk(isTRUE(pcv$exact) && all(dim(pcv$V) == c(4, 4)),
    "the pair covariance is the exact 4 x 4, not a block diagonal",
    "the pair covariance fell back to the block-diagonal form")
# WHICH pairs covary is not incidental -- it is the design. Group is BETWEEN
# participants, so "a x p" and "b x p" are estimated from disjoint participants
# and their cell coefficients are exactly independent. Condition is WITHIN, so
# "a x p" and "a x q" share every participant and must covary. An exact linear
# map reproduces both facts; a block-diagonal approximation gets the first one
# right by luck and the second one wrong.
cx <- function(c1, c2) {
  P <- e$dance_traj_pair_cov(co, match(c1, co$cells), match(c2, co$cells))
  max(abs(P$V[1:2, 3:4]))
}
chk(cx("a x p", "b x p") < 1e-12 && cx("a x p", "b x q") < 1e-12,
    "cells separated by the BETWEEN factor are exactly independent, as disjoint participants must be",
    "the map gives a non-zero covariance between cells with no participant in common")
chk(cx("a x p", "a x q") > 1e-6,
    sprintf("cells separated by the WITHIN factor covary (largest entry %.4f) -- assuming zero there would be wrong",
            cx("a x p", "a x q")),
    "the map gives zero covariance between two cells measured on the same participants")
# and the direction of the error, stated rather than assumed: Var(d) = V11 + V22 - 2C
pw2 <- e$dance_traj_pair_cov(co, match("a x p", co$cells), match("a x q", co$cells))
gd <- c(-1, 0, 1, 0)
v_exact <- as.numeric(t(gd) %*% pw2$V %*% gd)
Vbd <- pw2$V; Vbd[1:2, 3:4] <- 0; Vbd[3:4, 1:2] <- 0
v_drop <- as.numeric(t(gd) %*% Vbd %*% gd)
chk(abs(v_exact - v_drop) > 1e-8,
    sprintf("dropping it changes that contrast's variance: %.5f exact vs %.5f dropped, i.e. %s",
            v_exact, v_drop, if (v_drop < v_exact) "TOO NARROW" else "too wide"),
    "dropping the cross-cell block changed nothing, so this test proves nothing")

# ---------------------------------------------------------------- contrasts
cat("\n-- pairwise contrasts over the whole factorial ---------------------\n")
ct <- e$dance_traj_contrasts(fit, "level")
chk(isTRUE(ct$ok) && nrow(ct$table) == 6 && ct$n_cells == 4,
    "4 cells give all 6 pairs, not just the ones a 2-group design would have",
    sprintf("the contrast set is wrong: %d rows", nrow(ct$table %||% data.frame())))
chk(all(ct$table$p_adj >= ct$table$p_raw - 1e-12) && ct$adjust == "holm",
    "Holm adjustment is on by default and never lowers a p-value",
    "the multiplicity adjustment is missing or lowers a p-value")
ctn <- e$dance_traj_contrasts(fit, "level", adjust = "none")
chk(grepl("NO multiplicity adjustment", ctn$note) &&
      all(abs(ctn$table$p_adj - ctn$table$p_raw) < 1e-12),
    "and turning it off says plainly that the family is then wrong",
    "unadjusted contrasts came back without saying so")
# the level contrast is exactly linear, so it can be checked against the fit
i1 <- match("a x p", g$.cell); i2 <- match("b x p", g$.cell)
d1 <- e$dance_traj_design_rows(fit, g[i1, , drop = FALSE], fit$spec$t0) -
      e$dance_traj_design_rows(fit, g[i2, , drop = FALSE], fit$spec$t0)
man <- as.numeric(d1[, kp, drop = FALSE] %*% bv$beta[kp])
row <- ct$table[ct$table$cell1 == "a x p" & ct$table$cell2 == "b x p", ]
chk(nrow(row) == 1 && abs(row$estimate - man) < 1e-10,
    "a level contrast equals the linear combination it is defined as",
    "a level contrast does not match its own definition")

ca <- e$dance_traj_contrasts(fit, "amplitude")
chk(isTRUE(ca$ok) && nrow(ca$table) == 6 && all(is.finite(ca$table$se)),
    "amplitude contrasts come back for every pair with a finite SE",
    "amplitude contrasts failed")
cp <- e$dance_traj_contrasts(fit, "phase")
chk(isTRUE(cp$ok) && all(abs(cp$table$estimate) <= 12 + 1e-9),
    sprintf("phase contrasts are wrapped onto the effective period (max |diff| %.2f h <= 12)",
            max(abs(cp$table$estimate))),
    "a phase contrast left the half-period it must wrap onto")
chk(grepl("effective period", cp$unit) && grepl("amplitude difference", ca$unit),
    "each contrast set names the units its estimate is in",
    "a contrast set came back without units")

# Bingham again, this time on the contrast: no rhythm, so no phase contrast
fit_flat <- e$dance_traj_fit(e$dance_traj_spec(
  build(make(13, amp = 0, noise = 3)), 24, 1, "none"))
cpf <- e$dance_traj_contrasts(fit_flat, "phase")
chk(isTRUE(cpf$ok) && any(!cpf$table$defined) &&
      all(is.na(cpf$table$lo[!cpf$table$defined])),
    sprintf("with no rhythm, %d of %d phase contrasts are UNDEFINED and carry no interval",
            sum(!cpf$table$defined), nrow(cpf$table)),
    "a phase contrast was reported for cells whose amplitude may be zero")
chk(grepl("UNDEFINED", cpf$note %||% "") &&
      all(is.na(cpf$table$p_adj[!cpf$table$defined])),
    "and they are excluded from the multiplicity adjustment rather than counted as tests",
    "undefined contrasts were carried into the adjustment")

# ---------------------------------------------------------------- the adapter
cat("\n-- the cell-keyed adapter the existing plots consume ---------------\n")
gf <- e$dance_traj_group_fits(fit)
chk(length(gf) == 4 && identical(sort(names(gf)), sort(g$.cell)),
    "one entry per design cell, keyed by the cell label",
    "the adapter did not produce one entry per cell")
need <- c("group", "n", "mean_mesor", "intercept", "rhythm_adjusted_mean",
          "trend_coefs", "mean_coefs", "mean_amplitudes", "mean_acrophases_rad",
          "mean_acrophases_time", "sd_mesor", "sd_amplitudes")
chk(all(vapply(need, function(f) f %in% names(gf[[1]]), logical(1))),
    "every field the existing plot code reads is present",
    sprintf("missing: %s", paste(setdiff(need, names(gf[[1]])), collapse = ", ")))
chk(sum(vapply(gf, function(x) x$n, integer(1))) == 32 &&
      all(vapply(gf, function(x) x$n, integer(1)) == 8),
    "n is participants per cell, and the cells account for all 32 curves",
    "the per-cell participant counts do not reconcile")

# the coefficient vector must reconstruct the cell's own fitted curve
cellk <- "b x q"
cf <- gf[[cellk]]$mean_coefs
tt <- seq(0, 22, length.out = 40)
recon <- cf[1] + cf[2] * cos(2 * pi * tt / 24) + cf[3] * sin(2 * pi * tt / 24)
pred <- e$dance_traj_predict(fit, times = tt)$table
pk <- pred$fit[pred$cell == cellk]
chk(max(abs(recon - pk)) < 1e-8,
    sprintf("mean_coefs reconstructs that cell's fitted trajectory (max diff %.2e)",
            max(abs(recon - pk))),
    "the adapter's coefficient vector does not reproduce the model's own curve")

# THE RENAME THAT MUST NOT HAPPEN
chk(!is.null(gf[[1]]$se_amplitudes) && !is.null(gf[[1]]$sd_amplitudes) &&
      !isTRUE(all.equal(gf[[1]]$se_amplitudes, gf[[1]]$sd_amplitudes)),
    "precision (se_amplitudes) and dispersion (sd_amplitudes) are separate numbers",
    "the adapter reports the same number as both a standard error and a standard deviation")
chk(grepl("random-effects covariance", gf[[1]]$dispersion_source %||% ""),
    "and the dispersion says where it came from",
    "the dispersion came back unattributed")
chk(is.null(gf[[1]]$amp_arithmetic) && is.null(gf[[1]]$resultants) &&
      is.null(gf[[1]]$variance_decomp),
    "the three two-stage-only quantities are absent, not approximated",
    "the adapter invented a two-stage quantity a one-stage fit cannot supply")
gaps <- attr(gf, "unavailable")
chk(length(gaps) == 3 && all(nzchar(unlist(gaps))) &&
      grepl("E\\[v\\]", gaps$amp_arithmetic),
    "and each absence carries the reason it is absent",
    "the unavailable list is missing or empty")
chk(grepl("mixed-effects trajectory fit", attr(gf, "provenance") %||% ""),
    "the structure says which fit produced it, so it cannot be mistaken for the two-stage one",
    "the adapter's provenance is missing")

# a trend must survive the round trip too
fit_tr <- e$dance_traj_fit(e$dance_traj_spec(
  build(make(14, slope = 0.15)), 24, 1, "linear"))
gtr <- e$dance_traj_group_fits(fit_tr)
cf2 <- gtr[["a x p"]]$mean_coefs
recon2 <- cf2[1] + cf2[2] * (tt - fit_tr$spec$t0) +
          cf2[3] * cos(2 * pi * tt / 24) + cf2[4] * sin(2 * pi * tt / 24)
pred2 <- e$dance_traj_predict(fit_tr, times = tt)$table
chk(max(abs(recon2 - pred2$fit[pred2$cell == "a x p"])) < 1e-8,
    "with a linear trend the coefficient vector still reconstructs the curve, trend included",
    "the adapter's trend coefficient does not reproduce the model's curve")
chk(abs(gtr[["a x p"]]$rhythm_adjusted_mean - gtr[["a x p"]]$intercept) > 1e-6,
    "and the rhythm-adjusted mean is not the intercept when a trend is present",
    "the MESOR and the intercept came back identical under a trend")
}))

cat(sprintf("\n%s  (%d passed, %d failed)\n",
            if (fail == 0) "Trajectory band/contrast tests PASSED" else "FAILURES",
            pass, fail))
if (fail > 0) quit(status = 1)
