# ==============================================================================
# tests/testthat/test-apa-report.R — the publication report
#
# The report's job is to state, in prose a journal will accept, both the
# analytic choices and the results. Two things can go wrong with such a
# document, and neither is a crash:
#
#   * it describes an analysis that was not run, or omits one that was;
#   * it formats a statistic in a way APA does not use, or -- worse -- states
#     a descriptive quantity as though it were inferential.
#
# The end-to-end generation is exercised in tests/reactive_smoke_test.R, which
# builds a real session and reads the document it produces. These are the unit
# checks on the formatting rules and on the promises the text makes.
# ==============================================================================
.libPaths(c("~/Rlib", .libPaths()))
app_dir <- if (dir.exists("server")) "." else if (dir.exists("../../server")) "../.." else "FCK"
source(file.path(app_dir, "server/93_apa_report.R"), local = TRUE)

# --------------------------------------------------------------- APA 7 rules -
test_that("a quantity that cannot exceed 1 loses its leading zero", {
  # APA 7 sec. 6.36
  expect_identical(fck_apa_num(0.427, 3, bounded = TRUE), ".427")
  expect_identical(fck_apa_num(-0.009, 3, bounded = TRUE), "-.009")
  expect_identical(fck_apa_num(0.427, 2, bounded = FALSE), "0.43")
  expect_identical(fck_apa_num(12.3456, 2), "12.35")
  expect_identical(fck_apa_num(NA), "--")
  expect_identical(fck_apa_num(NULL), "--")
  expect_identical(fck_apa_num(Inf), "--")
  # vectorised, because variance tables are built from one call
  expect_identical(fck_apa_num(c(0.1, 0.25), 2, TRUE), c(".10", ".25"))
})

test_that("p-values follow APA 7 sec. 6.42", {
  expect_identical(fck_apa_p(0.0432), "*p* = .043")
  expect_identical(fck_apa_p(0.0004), "*p* < .001")
  expect_identical(fck_apa_p(0.9999), "*p* > .999")
  expect_identical(fck_apa_pval(0.0432), ".043")
  expect_identical(fck_apa_pval(0.0004), "< .001")
  expect_identical(fck_apa_pval(NA), "--")
  # never a leading zero, at any precision
  for (p in c(.5, .05, .005, .0005))
    expect_false(grepl("^0\\.", fck_apa_pval(p, 4)))
})

test_that("F is reported with exact degrees of freedom", {
  expect_identical(fck_apa_F(12.876, 3, 122), "*F*(3, 122) = 12.88")
  expect_identical(fck_apa_F(12.876, 3, 122, 0.0004),
                   "*F*(3, 122) = 12.88, *p* < .001")
  expect_identical(fck_apa_F(NA, 3, 122), "*F* undefined")
})

test_that("a bounded RANGE drops the leading zero at both ends", {
  expect_identical(fck_apa_range(c(0.007, 0.241), 3, bounded = TRUE), ".007 to .241")
  expect_identical(fck_apa_range(c(0.007, 0.241), 3), "0.007 to 0.241")
  expect_identical(fck_apa_range(c(NA, NaN)), "--")
})

test_that("eponymous corrections are written as proper names", {
  expect_identical(fck_apa_method("bonferroni"), "Bonferroni")
  expect_identical(fck_apa_method("holm"), "Holm")
  expect_identical(fck_apa_method("BH"), "Benjamini-Hochberg")
  expect_identical(fck_apa_method("fdr"), "Benjamini-Hochberg")
  expect_identical(fck_apa_method("none"), "no")
  expect_identical(fck_apa_method("something-else"), "something-else")
})

test_that("M and SD are reported together, and drop empty input", {
  expect_identical(fck_apa_M(c(1, 2, 3)), "*M* = 2.00, *SD* = 1.00")
  expect_identical(fck_apa_M(c(NA, NaN)), "--")
})

# ------------------------------------------------------------ the structure --
test_that("a report on an empty session says so and stops", {
  v <- list(data = NULL)
  md <- fck_apa_report(v, list())
  txt <- paste(md, collapse = "\n")
  expect_true(grepl("No data had been imported", txt, fixed = TRUE))
  # and it does not invent Methods or Results for analyses that cannot exist.
  # (The how-to-use note legitimately contains the WORD "Results"; what must be
  # absent is the SECTION.)
  expect_false(grepl("Functional ANOVA", txt, fixed = TRUE))
  expect_false(any(grepl("^## Results$", md)))
  expect_false(any(grepl("^## Statistical analysis$", md)))
})

test_that("a report on data with no analyses describes the data and says so", {
  set.seed(1)
  v <- list(data = matrix(rnorm(20 * 6), 20, 6),
            time_labels = sprintf("%02d:00", 0:5),
            subject_ids = paste0("S", 1:20))
  md <- fck_apa_report(v, list())
  txt <- paste(md, collapse = "\n")
  expect_true(grepl("20 curves observed at 6 time points", txt, fixed = TRUE))
  expect_true(grepl("00:00 to 05:00", txt, fixed = TRUE))
  expect_true(grepl("No analysis had been run", txt, fixed = TRUE))
  expect_true(grepl("Reproducibility", txt, fixed = TRUE))
})

test_that("repeated observations of one participant are flagged, not hidden", {
  set.seed(2)
  v <- list(data = matrix(rnorm(10 * 4), 10, 4),
            subject_ids = rep(paste0("S", 1:5), 2))
  txt <- paste(fck_apa_report(v, list()), collapse = "\n")
  expect_true(grepl("5 distinct participants", txt, fixed = TRUE))
  expect_true(grepl("anticonservative", txt, fixed = TRUE))
})

test_that("the permutation resolution floor is stated with the result", {
  set.seed(3)
  v <- list(data = matrix(rnorm(12 * 5), 12, 5),
            fanova_results = list(
              design = "between", n_groups = 2, groups = c("A", "B"),
              group_sizes = c(6L, 6L), n_permutations = 100L, alpha = 0.05,
              F_stat = rep(2, 5), p_values_pointwise = rep(.02, 5),
              p_values_adjusted = rep(.02, 5), sig_regions = rep(TRUE, 5),
              eta_squared = rep(.2, 5), eta_squared_type = "classical",
              df_between = 1, df_within = 10, time_points = seq(0, 1, length.out = 5),
              L2_stat = 1, p_value_L2 = .0099, n_undefined = 0))
  txt <- paste(fck_apa_report(v, list()), collapse = "\n")
  # 1/(100+1) = .0099
  expect_true(grepl("smallest attainable value here is .0099", txt, fixed = TRUE))
  expect_true(grepl("*B* = 100", txt, fixed = TRUE))
})

test_that("every Results subsection states what its numbers do not establish", {
  set.seed(4)
  v <- list(data = matrix(rnorm(12 * 5), 12, 5),
            smooth_fit_metrics = list(method = "auto", n_basis = 8,
                                      lambda = 1e-4, mean_r_squared = .97,
                                      mean_rmse = .11, mean_df = 6,
                                      time_axis = "column index"),
            clustering_results = list(k = 3L, cluster_sizes = c(4L, 4L, 4L),
                                      r_squared = .6, silhouette_width = .41,
                                      ch_index = 12.2, method = "kmeans",
                                      method_label = "functional k-means",
                                      standardize = "none"))
  md  <- fck_apa_report(v, list())
  txt <- paste(md, collapse = "\n")
  n_h3 <- sum(grepl("^### ", md))
  # count the caveat PARAGRAPHS, not the phrase -- the how-to-use note names it
  n_caveat <- sum(grepl("^\\*What these numbers do not establish\\.\\*", md))
  expect_gt(n_h3, 0)
  expect_equal(n_caveat, n_h3)
  # the specific claims this audit had to remove must not come back
  expect_true(grepl("k-means partitions any data set", txt, fixed = TRUE))
  expect_false(grepl("variance explained by phase", txt, fixed = TRUE) &&
               !grepl("Do not report", txt, fixed = TRUE))
})

test_that("registration diagnostics refuse the two claims P5.2 and P5.3 removed", {
  set.seed(5)
  v <- list(data = matrix(rnorm(12 * 5), 12, 5),
            warping_results = list(
              method = "linear_shift", boundary = "periodic wrap",
              fit_statistics = list(summary = list(
                dispersion_reduction = .42, total_dispersion_pre = 10,
                total_dispersion_post = 5.8, mean_phase_displacement = .08))))
  txt <- paste(fck_apa_report(v, list()), collapse = "\n")
  expect_true(grepl("not** an", txt, fixed = TRUE))
  expect_true(grepl("amplitude/phase variance decomposition", txt, fixed = TRUE))
  expect_true(grepl("no likelihood to penalise", txt, fixed = TRUE))
  expect_false(grepl("AIC", txt, fixed = TRUE))
  expect_false(grepl("BIC", txt, fixed = TRUE))
})

test_that("the amplitude-leakage signature is flagged in the report too", {
  v <- list(data = matrix(1, 4, 4),
            warping_results = list(method = "parametric", family = "logistic",
              fit_statistics = list(summary = list(
                dispersion_reduction = .28, total_dispersion_pre = 10,
                total_dispersion_post = 7.2, mean_phase_displacement = .014))))
  txt <- paste(fck_apa_report(v, list()), collapse = "\n")
  expect_true(grepl("Caution", txt, fixed = TRUE))
  expect_true(grepl("absorbed *amplitude* differences", txt, fixed = TRUE))
})

# ---------------------------------------------------------------- HTML ------
test_that("the HTML rendering keeps tables and resolves escaped asterisks", {
  md <- c("# Title", "", "*Generated by F\\*CK on 2026-01-01.*", "",
          "> A note with **bold**.", "",
          "| A | B |", "| :--- | ---: |", "| x | 1 |", "", "---", "",
          "A paragraph with *emphasis* and `code`.")
  h <- paste(fck_apa_html(md, "T"), collapse = "\n")
  expect_true(grepl("<!DOCTYPE html>", h, fixed = TRUE))
  expect_true(grepl("<h1>Title</h1>", h, fixed = TRUE))
  expect_true(grepl("<table>", h, fixed = TRUE))
  expect_true(grepl("<blockquote>", h, fixed = TRUE))
  expect_true(grepl("<strong>bold</strong>", h, fixed = TRUE))
  expect_true(grepl("<em>emphasis</em>", h, fixed = TRUE))
  expect_true(grepl("<code>code</code>", h, fixed = TRUE))
  expect_true(grepl("<hr>", h, fixed = TRUE))
  # the escaped asterisk survives as a literal and does not break the italics
  expect_true(grepl("F*CK", h, fixed = TRUE))
  expect_false(grepl("F\\CK", h, fixed = TRUE))
  # no Markdown table syntax leaks through
  expect_false(grepl("| :---", h, fixed = TRUE))
  expect_false(grepl("---: |", h, fixed = TRUE))
})

test_that("the HTML escapes angle brackets but keeps the tags the report uses", {
  h <- paste(fck_apa_html(c("*R*<sup>2</sup> and a <script>alert(1)</script>")),
             collapse = "\n")
  expect_true(grepl("<sup>2</sup>", h, fixed = TRUE))
  expect_true(grepl("&lt;script&gt;", h, fixed = TRUE))
  expect_false(grepl("<script>", h, fixed = TRUE))
})
