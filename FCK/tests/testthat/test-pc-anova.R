# ==============================================================================
# tests/testthat/test-pc-anova.R — per-component group tests on fPCA scores
#
# The things that are easy to get wrong here and that these pin down:
#   * TWO multiplicity families. Correcting the pairwise tests and forgetting
#     the k omnibus tests is the usual error; the across-component correction
#     must actually change the omnibus p-values.
#   * The post-hoc gate must run on the ACROSS-COMPONENT-ADJUSTED p, so a
#     component that only survives before that correction gets no pairwise
#     table lending it credibility.
#   * Welch must be selected when the variances differ, not assumed away.
#   * Tukey and Games-Howell are family-wise BY CONSTRUCTION, so p.adjust must
#     not be applied to them a second time.
#   * Eigenvalue separation must flag near-ties, where "a difference on PC2" is
#     a difference on an arbitrary rotation of PC2 and PC3.
#   * Repeated curves must be detected, because a between-groups test on them
#     is anticonservative and nothing else in the pipeline would notice.
# ==============================================================================

`%||%` <- function(a, b) if (is.null(a)) b else a
app_dir <- if (dir.exists("server")) "." else if (dir.exists("../../server")) "../.." else "FCK"
source(file.path(app_dir, "server/08_helpers_cosinor.R"))
source(file.path(app_dir, "server/09_helpers_pcanova.R"))

mk <- function(n = c(60, 60, 60), shift = c(0, 0, 0), sd = c(1, 1, 1), seed = 1) {
  set.seed(seed)
  g <- factor(rep(c("A", "B", "C")[seq_along(n)], n))
  x <- unlist(Map(function(k, m, s) rnorm(k, m, s), n, shift, sd))
  list(x = x, g = g)
}

# ------------------------------------------------------------------- omnibus
test_that("the omnibus returns Fisher, Welch and Kruskal-Wallis together", {
  d <- mk(shift = c(0, 0.6, 1.2))
  o <- fck_pc_omnibus(d$x, d$g)
  expect_equal(o$k, 3)
  expect_true(all(is.finite(c(o$F, o$p, o$welch_F, o$welch_p, o$kw_chisq, o$kw_p))))
  # against R's own implementations
  expect_equal(o$p, stats::anova(stats::lm(d$x ~ d$g))$`Pr(>F)`[1], tolerance = 1e-10)
  expect_equal(o$welch_p, stats::oneway.test(d$x ~ d$g, var.equal = FALSE)$p.value,
               tolerance = 1e-10)
})

test_that("omega2 is never larger than eta2", {
  for (s in 1:20) {
    d <- mk(shift = c(0, runif(1, 0, 1), runif(1, 0, 1)), seed = s)
    o <- fck_pc_omnibus(d$x, d$g)
    expect_lte(o$omega2, o$eta2)
  }
})

test_that("a null effect gives a large p and an omega2 near (or below) zero", {
  d <- mk(shift = c(0, 0, 0), seed = 1)
  o <- fck_pc_omnibus(d$x, d$g)
  expect_gt(o$p, 0.30)
  expect_lt(o$omega2, 0.02)
  # omega2 is legitimately NEGATIVE under a true null -- it subtracts the df
  # penalty from the between-groups SS -- and must not be clamped to 0, because
  # clamping would make "no effect" indistinguishable from "a trivial effect".
  expect_lt(o$omega2, o$eta2)
  expect_gte(o$eta2, 0)
})

# --------------------------------------------------------- assumption checks
test_that("Brown-Forsythe detects unequal variances and misses equal ones", {
  het <- mk(sd = c(1, 1, 4), seed = 3)
  hom <- mk(sd = c(1, 1, 1), seed = 1)
  expect_lt(fck_brown_forsythe(het$x, het$g)$p, 0.01)
  expect_gt(fck_brown_forsythe(hom$x, hom$g)$p, 0.05)
})

test_that("normality reporting separates severity from Shapiro's p", {
  set.seed(5)
  g <- factor(rep(c("A", "B"), each = 400))
  # mildly non-normal: Shapiro will reject at n = 800, but ANOVA is fine
  mild <- fck_normality_check(rnorm(800)^2 * sign(rnorm(800)), g)
  expect_true(is.finite(mild$skew))
  # severely non-normal
  sev <- fck_normality_check(c(rnorm(790), rnorm(10, 0, 40)), g)
  expect_true(sev$severe)
})

# -------------------------------------------------------- the Welch decision
test_that("Welch is used automatically when the variances differ", {
  d <- mk(n = c(300, 60), shift = c(0, 0.5), sd = c(1, 4), seed = 7)
  S <- matrix(d$x, ncol = 1)
  r <- fck_pc_anova_all(S, d$g, eigenvalues = c(1), across_pc_correction = "none")
  o <- r$omnibus[[1]]
  expect_equal(r$p_omnibus[1], o$welch_p, tolerance = 1e-12)
  # and Fisher is NOT what was carried forward
  expect_false(isTRUE(all.equal(r$p_omnibus[1], o$p, tolerance = 1e-6)))
})

test_that("the user can force Fisher or Welch", {
  d <- mk(n = c(300, 60), shift = c(0, 0.5), sd = c(1, 4), seed = 7)
  S <- matrix(d$x, ncol = 1)
  rf <- fck_pc_anova_all(S, d$g, across_pc_correction = "none", use_welch = FALSE)
  rw <- fck_pc_anova_all(S, d$g, across_pc_correction = "none", use_welch = TRUE)
  expect_equal(rf$p_omnibus[1], rf$omnibus[[1]]$p, tolerance = 1e-12)
  expect_equal(rw$p_omnibus[1], rw$omnibus[[1]]$welch_p, tolerance = 1e-12)
})

# ------------------------------------------------------ the TWO families
test_that("the across-component correction actually adjusts the omnibus p's", {
  set.seed(21)
  g <- factor(rep(c("A", "B", "C"), each = 80))
  S <- matrix(rnorm(240 * 6), ncol = 6)          # pure noise, 6 components
  raw <- fck_pc_anova_all(S, g, across_pc_correction = "none")
  adj <- fck_pc_anova_all(S, g, across_pc_correction = "holm")
  bon <- fck_pc_anova_all(S, g, across_pc_correction = "bonferroni")
  expect_equal(raw$p_omnibus_adj, raw$p_omnibus, tolerance = 1e-12)
  expect_true(all(adj$p_omnibus_adj >= adj$p_omnibus - 1e-12))
  expect_true(all(bon$p_omnibus_adj >= adj$p_omnibus_adj - 1e-12))  # Bonferroni >= Holm
  expect_equal(bon$p_omnibus_adj, pmin(1, 6 * bon$p_omnibus), tolerance = 1e-10)
})

test_that("the post-hoc gate runs on the ADJUSTED omnibus p, not the raw one", {
  set.seed(1)
  g <- factor(rep(c("A", "B"), each = 90))
  S <- matrix(rnorm(180 * 8), ncol = 8)
  # An effect deliberately sized to land BETWEEN the two thresholds: raw
  # p = 0.028 (passes a 0.05 gate) but Bonferroni over 8 components gives
  # p = 0.22 (fails it). This is exactly the case the gate exists to catch.
  S[g == "B", 1] <- S[g == "B", 1] + 0.41
  raw <- fck_pc_anova_all(S, g, across_pc_correction = "none", posthoc_gate = 0.05)
  bon <- fck_pc_anova_all(S, g, across_pc_correction = "bonferroni", posthoc_gate = 0.05)
  expect_lt(raw$p_omnibus[1], 0.05)
  expect_gt(raw$p_omnibus[1], 0.01)
  expect_gt(bon$p_omnibus_adj[1], 0.05)
  expect_false(is.null(raw$posthoc[[1]]))   # gate opened on the raw p
  expect_null(bon$posthoc[[1]])             # and closed on the adjusted one
})

test_that("components above the gate get no post-hoc table at all", {
  set.seed(8)
  g <- factor(rep(c("A", "B", "C"), each = 70))
  S <- matrix(rnorm(210 * 4), ncol = 4)
  S[g == "C", 1] <- S[g == "C", 1] + 2       # PC1 real, PC2-4 noise
  r <- fck_pc_anova_all(S, g, across_pc_correction = "holm", posthoc_gate = 0.05)
  expect_false(is.null(r$posthoc[[1]]))
  expect_true(all(vapply(r$posthoc[2:4], is.null, logical(1))))
})

# ------------------------------------------------------------- the post-hocs
test_that("every correction is monotone in the raw p and never smaller", {
  d <- mk(n = c(70, 70, 70, 70), shift = c(0, 0.4, 0.8, 1.2), sd = c(1, 1, 1, 1), seed = 12)
  d$g <- factor(rep(c("A", "B", "C", "D"), each = 70))
  for (m in c("bonferroni", "holm", "hochberg", "hommel", "BH", "BY")) {
    ph <- fck_pc_posthoc(d$x, d$g, method = m)
    expect_true(all(ph$p_adj >= ph$p_raw - 1e-12), info = m)
    expect_true(all(ph$p_adj <= 1 + 1e-12), info = m)
  }
  none <- fck_pc_posthoc(d$x, d$g, method = "none")
  expect_equal(none$p_adj, none$p_raw, tolerance = 1e-12)
})

test_that("Tukey and Games-Howell are not double-adjusted", {
  d <- mk(n = c(80, 80, 80), shift = c(0, 0.5, 1), seed = 15)
  for (m in c("tukey", "games-howell")) {
    ph <- fck_pc_posthoc(d$x, d$g, method = m)
    expect_true(all(ph$already_adjusted))
    expect_equal(ph$p_adj, ph$p_raw, tolerance = 1e-12)   # p.adjust NOT applied again
  }
})

test_that("Tukey reproduces TukeyHSD on a balanced equal-variance design", {
  d <- mk(n = c(60, 60, 60), shift = c(0, 0.5, 1.1), sd = c(1, 1, 1), seed = 31)
  ph <- fck_pc_posthoc(d$x, d$g, method = "tukey")
  ref <- TukeyHSD(aov(d$x ~ d$g))[[1]]
  # rows are "B-A","C-A","C-B" in TukeyHSD; ours are A-B, A-C, B-C (sign flipped)
  key <- paste0(ph$b, "-", ph$a)
  expect_equal(sort(ph$p_raw), sort(as.numeric(ref[, "p adj"])), tolerance = 1e-6)
  for (i in seq_len(nrow(ph)))
    expect_equal(-ph$diff[i], unname(ref[key[i], "diff"]), tolerance = 1e-8)
})

test_that("Games-Howell differs from Tukey when the variances differ", {
  d <- mk(n = c(200, 40, 40), shift = c(0, 0.5, 1), sd = c(1, 5, 5), seed = 17)
  tk <- fck_pc_posthoc(d$x, d$g, method = "tukey")
  gh <- fck_pc_posthoc(d$x, d$g, method = "games-howell")
  expect_false(isTRUE(all.equal(tk$p_raw, gh$p_raw, tolerance = 1e-3)))
  # Games-Howell's intervals are wider where the small, noisy groups are involved
  expect_gt(mean(gh$ci_hi - gh$ci_lo), mean(tk$ci_hi - tk$ci_lo))
})

test_that("every contrast carries an interval and an effect size", {
  d <- mk(n = c(50, 50, 50), shift = c(0, 1, 2), seed = 23)
  ph <- fck_pc_posthoc(d$x, d$g, method = "holm")
  expect_equal(nrow(ph), 3)                      # 3 choose 2
  expect_true(all(is.finite(ph$hedges_g)))
  expect_true(all(ph$ci_lo < ph$ci_hi))
  # the interval must bracket the difference
  expect_true(all(ph$ci_lo <= ph$diff & ph$diff <= ph$ci_hi))
  # A vs C is the biggest gap and should have the biggest |g|
  expect_equal(which.max(abs(ph$hedges_g)),
               which(ph$a == "A" & ph$b == "C"))
})

test_that("the number of pairwise tests is k choose 2", {
  for (k in 2:5) {
    n <- rep(40, k)
    g <- factor(rep(LETTERS[1:k], n))
    set.seed(k); x <- rnorm(sum(n))
    expect_equal(nrow(fck_pc_posthoc(x, g, "holm")), choose(k, 2))
  }
})

# ------------------------------------------------------- eigenvalue separation
test_that("near-tied eigenvalues are flagged as poorly separated", {
  sp <- fck_pc_separation(c(10, 5, 4.9, 1))
  expect_equal(sp$stability[1], "well separated")
  expect_true(grepl("POORLY", sp$stability[2]))
  expect_true(grepl("POORLY", sp$stability[3]))
})

test_that("a clearly ordered spectrum is not flagged", {
  sp <- fck_pc_separation(c(100, 20, 5, 1))
  expect_true(all(sp$stability == "well separated"))
})

# --------------------------------------------------------- independence check
test_that("repeated curves are detected", {
  ic <- fck_pc_independence_check(c(1, 1, 2, 3, 3, 3))
  expect_true(ic$checked)
  expect_false(ic$ok)
  expect_equal(ic$n_rows, 6)
  expect_equal(ic$n_subjects, 3)
  expect_equal(ic$n_repeated, 3)
  expect_equal(ic$max_per_subject, 3)
})

test_that("one row per participant passes, and no id means no claim", {
  expect_true(fck_pc_independence_check(1:10)$ok)
  expect_false(fck_pc_independence_check(NULL)$checked)
})

# ------------------------------------------------------------- the whole table
test_that("the full runner finds a planted effect and rejects the noise", {
  set.seed(2026)
  g <- factor(rep(c("YOUTH", "ADULT", "ELDERLY"), c(300, 200, 60)))
  S <- matrix(rnorm(560 * 5), ncol = 5)
  S[g == "ELDERLY", 2] <- S[g == "ELDERLY", 2] + 1.5    # a real effect on PC2 only
  r <- fck_pc_anova_all(S, g, eigenvalues = c(40, 20, 8, 3, 1),
                        varprop = c(.5, .25, .1, .1, .05),
                        across_pc_correction = "holm",
                        posthoc_method = "games-howell")
  expect_equal(r$k, 5)
  expect_lt(r$p_omnibus_adj[2], 0.001)
  expect_true(all(r$p_omnibus_adj[c(1, 3, 4, 5)] > 0.05))
  expect_false(is.null(r$posthoc[[2]]))
  expect_true(all(vapply(r$posthoc[c(1, 3, 4, 5)], is.null, logical(1))))
  # the ELDERLY contrasts are the significant ones
  ph <- r$posthoc[[2]]
  eld <- ph$a == "ELDERLY" | ph$b == "ELDERLY"
  expect_true(all(ph$p_adj[eld] < 0.05))
  expect_true(all(ph$p_adj[!eld] > 0.05))
})

test_that("an untestable component keeps its slot instead of renumbering the rest", {
  # `l[[j]] <- NULL` DELETES element j in R. Building the per-component lists
  # that way meant one untestable component shifted every later component down
  # by one, so PC4's result would be reported as PC3's. This pins the fix.
  set.seed(2)
  g <- factor(rep(c("A", "B"), each = 20))
  S <- cbind(rnorm(40), rep(NA_real_, 40), rnorm(40))
  r <- fck_pc_anova_all(S, g, across_pc_correction = "none")
  expect_false(is.null(r))
  expect_equal(r$k, 3)
  expect_equal(length(r$omnibus), 3)     # not 2
  expect_equal(length(r$variance), 3)
  expect_equal(length(r$posthoc), 3)
  expect_false(is.null(r$omnibus[[1]]))
  expect_null(r$omnibus[[2]])            # the slot survives, holding NULL
  expect_false(is.null(r$omnibus[[3]]))  # and PC3 is still PC3
})

test_that("a single-level grouping variable yields nothing rather than nonsense", {
  g <- factor(rep("A", 30))
  expect_null(fck_pc_omnibus(rnorm(30), g))
})
