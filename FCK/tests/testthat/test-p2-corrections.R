# ==============================================================================
# tests/testthat/test-p2-corrections.R
#
# P2: engineering. None of this changes a number, so the tests are mostly about
# preventing drift -- except the vectorisation, which MUST be arithmetically
# identical to the loop it replaced and is checked against it directly.
# ==============================================================================
.libPaths(c("~/Rlib", .libPaths()))
`%||%` <- function(a, b) if (is.null(a)) b else a
app_dir <- if (dir.exists("server")) "." else if (dir.exists("../../server")) "../.." else "FCK"
source(file.path(app_dir, "server/01b_kernel_contracts.R"))

# Strip comments so the guards read CODE, not prose -- the fixes deliberately
# quote the removed lines in their comments. A naive sub("#.*$") also eats code
# that follows a "#" INSIDE a string literal (e.g. add("#  ", R.version.string)),
# which would make a forbidden-string guard pass for the wrong reason. Only strip
# a "#" that is not inside quotes.
code_of <- function(f) {
  ln <- readLines(file.path(app_dir, f), warn = FALSE)
  strip <- function(l) {
    ch <- strsplit(l, "", fixed = TRUE)[[1]]
    if (!length(ch)) return(l)
    inq <- ""; out <- character(0)
    for (i in seq_along(ch)) {
      c1 <- ch[i]
      esc <- i > 1 && ch[i - 1] == "\\"
      if (nzchar(inq)) {
        if (c1 == inq && !esc) inq <- ""
      } else if (c1 %in% c("\"", "'")) {
        inq <- c1
      } else if (c1 == "#") {
        break
      }
      out <- c(out, c1)
    }
    paste(out, collapse = "")
  }
  paste(vapply(ln, strip, character(1), USE.NAMES = FALSE), collapse = "\n")
}

# --------------------------------------------- P2.5 vectorised permutation
test_that("the vectorised permutation statistic equals the loop it replaced", {
  set.seed(3)
  n_time <- 40; n_curves <- 60; n_groups <- 3
  curves <- matrix(rnorm(n_time * n_curves), n_time, n_curves)
  group_labels <- factor(rep(c("a", "b", "c"), length.out = n_curves))
  groups <- levels(group_labels)
  overall_mean <- rowMeans(curves)
  set.seed(9); perm_labels <- sample(group_labels)

  # the loop that used to be in server/50_fanova.R
  pm_old <- matrix(NA, n_time, n_groups)
  for (i in 1:n_groups)
    pm_old[, i] <- rowMeans(curves[, which(perm_labels == groups[i]), drop = FALSE])
  SSB_old <- numeric(n_time); SSW_old <- numeric(n_time)
  for (t in 1:n_time) for (i in 1:n_groups) {
    pidx <- which(perm_labels == groups[i])
    SSB_old[t] <- SSB_old[t] + length(pidx) * (pm_old[t, i] - overall_mean[t])^2
    for (j in pidx) SSW_old[t] <- SSW_old[t] + (curves[t, j] - pm_old[t, i])^2
  }

  # the form that ships
  gidx <- match(as.character(perm_labels), groups)
  pm <- vapply(seq_len(n_groups), function(i)
    rowMeans(curves[, gidx == i, drop = FALSE]), numeric(n_time))
  npg <- tabulate(gidx, nbins = n_groups)
  SSB_new <- as.vector((pm - overall_mean)^2 %*% npg)
  SSW_new <- rowSums((curves - pm[, gidx, drop = FALSE])^2)

  expect_equal(SSB_old, SSB_new, tolerance = 1e-10)
  expect_equal(SSW_old, SSW_new, tolerance = 1e-10)
  expect_equal(pm_old, pm, tolerance = 1e-12)
})

test_that("the four-deep interpreted loop is gone from the source", {
  src <- code_of("server/50_fanova.R")
  expect_false(grepl("for(j in perm_idx) {", src, fixed = TRUE))
  expect_true(grepl("SSW_perm <- rowSums((curves - perm_means[, gidx, drop = FALSE])^2)",
                    src, fixed = TRUE))
})

# ------------------------------------------------------- P2.2 no self-install
test_that("the app does not install packages while starting", {
  src <- code_of("app.R")
  expect_false(grepl("install_if_missing", src, fixed = TRUE))
  # install.packages may appear only inside the message telling the user what to
  # run, never as a call the app makes. A top-level call starts a line.
  ln <- readLines(file.path(app_dir, "app.R"), warn = FALSE)
  expect_equal(sum(grepl("^\\s*install\\.packages\\(", ln)), 0L)
  expect_equal(sum(grepl("^\\s*renv::(install|restore)\\(", ln)), 0L)
  expect_true(grepl("missing_required", src, fixed = TRUE))
  expect_true(grepl("cannot start: required packages are missing", src, fixed = TRUE))
})

test_that("a lockfile bootstrap exists and refuses to install renv for you", {
  expect_true(file.exists(file.path(app_dir, "tools/renv_bootstrap.R")))
  src <- readLines(file.path(app_dir, "tools/renv_bootstrap.R"), warn = FALSE)
  expect_true(any(grepl("renv::snapshot", src, fixed = TRUE)))
  expect_true(any(grepl("deliberately does not install it for you", src, fixed = TRUE)))
})

# ------------------------------------- P2.1 exported code calls the real kernel
test_that("the fANOVA export deparses the app's own estimator", {
  src <- code_of("server/90_export.R")
  # the reimplementation that could disagree with the app
  expect_false(grepl("aov_result <- summary(aov(curves_eval[, t] ~ group_labels))",
                     src, fixed = TRUE))
  expect_false(grepl("rm_result <- rmfanova(formula = value ~ condition", src, fixed = TRUE))
  # The live function objects, written out. P3.5 replaced the three separate
  # deparse() calls with emit_kernel(), which deparses the SAME objects and
  # additionally emits their helper dependency closure -- so the guard now
  # checks that each kernel is emitted through it, and that deparse of a live
  # object is still what does the writing.
  expect_true(grepl('emit_kernel("perform_functional_anova")', src, fixed = TRUE))
  expect_true(grepl('emit_kernel("perform_rm_fanova")', src, fixed = TRUE))
  expect_true(grepl('emit_kernel("fit_cosinor")', src, fixed = TRUE))
  expect_true(grepl('emit_kernel("fck_fit_fosr_ols")', src, fixed = TRUE))
  expect_true(grepl("paste(deparse(obj), collapse", src, fixed = TRUE))
  # the closure is computed from the live objects, not from a list that can go
  # stale -- that is the whole point of P3.5
  expect_true(grepl("codetools::findGlobals(obj, merge = TRUE)", src, fixed = TRUE))
  expect_true(grepl("kernel_env  <- environment(generate_analysis_code)", src, fixed = TRUE))
})

test_that("the exported script records and checks the environment it ran in", {
  src <- code_of("server/90_export.R")
  expect_true(grepl("fck_recorded_versions", src, fixed = TRUE))
  expect_true(grepl("fck_check_env", src, fixed = TRUE))
  expect_true(grepl("R.version.string", src, fixed = TRUE))
  # the check reports and continues; it must not abort someone's rerun
  expect_false(grepl("stop('Package versions differ", src, fixed = TRUE))
})

# ------------------------------------------------------- P2.3 kernel contracts
test_that("every statistical kernel has a recorded contract", {
  needed <- c("fit_cosinor", "fit_cosinor_nonlinear", "fck_auto_lambda",
              "perform_functional_anova", "perform_rm_fanova",
              "linear_shift_alignment", "fit_fosr",
              "run_functional_clustering")
  expect_true(all(needed %in% names(FCK_CONTRACTS)))
  # orientation is stated for every one: the app mixes subjects x time and
  # time x subjects, which is the easiest silent transpose in the codebase
  for (k in names(FCK_CONTRACTS))
    expect_true("orientation" %in% names(FCK_CONTRACTS[[k]]) ||
                "estimand" %in% names(FCK_CONTRACTS[[k]]), info = k)
  expect_error(fck_contract("no_such_kernel"), "No contract recorded")
})

test_that("the contracts record the facts the audit established", {
  expect_true(grepl("14.6%", FCK_CONTRACTS$fit_cosinor_nonlinear$zero_amp))
  expect_true(grepl("22.8%", FCK_CONTRACTS$perform_rm_fanova$estimator))
  expect_true(grepl("does NOT call rmfanova", FCK_CONTRACTS$perform_rm_fanova$package))
  expect_true(grepl("p_value_L2 is NA", FCK_CONTRACTS$perform_rm_fanova$global))
  expect_true(grepl("does not search", FCK_CONTRACTS$fck_auto_lambda$warning))
  expect_true(grepl("Not simultaneous", FCK_CONTRACTS$perform_functional_anova$intervals, fixed = TRUE))
  expect_null(FCK_CONTRACTS$fit_sofr)   # the module was removed
})

# ------------------------------------------------- P2.6 rmfanova global test
test_that("the RM design is built as rmfanova documents it", {
  source(file.path(app_dir, "server/09b_helpers_rmfanova.R"), local = TRUE)
  set.seed(5); n <- 10; p <- 12; l <- 3
  curves <- matrix(rnorm(p * n * l), p, n * l)          # time x curves
  sid <- factor(rep(1:n, times = l))
  rmf <- factor(rep(c("A", "B", "C"), each = n))

  d <- fck_rm_design(curves, sid, rmf)
  expect_equal(length(d$x), l)                          # a list of l matrices
  expect_true(all(vapply(d$x, nrow, 1L) == n))          # each n x p
  expect_true(all(vapply(d$x, ncol, 1L) == p))
  # the SAME subject must occupy the same row in every condition
  expect_true(all(vapply(d$x, function(m) identical(rownames(m), rownames(d$x[[1]])),
                         logical(1))))
  expect_equal(d$n_complete, n)

  # a subject missing one condition is dropped, and named
  drop <- which(sid == 3 & rmf == "C")
  d2 <- fck_rm_design(curves[, -drop], sid[-drop], rmf[-drop])
  expect_equal(d2$n_complete, n - 1)
  expect_true("3" %in% d2$dropped)
  expect_false("3" %in% rownames(d2$x[[1]]))

  # too few complete subjects is a refusal with a reason, not a partial answer
  d3 <- fck_rm_design(curves[, 1:4], sid[1:4], rmf[1:4])
  expect_null(d3$x)
  expect_true(nzchar(d3$reason))
})

test_that("only the calibrated subset is shown by default", {
  source(file.path(app_dir, "server/09b_helpers_rmfanova.R"), local = TRUE)
  expect_equal(FCK_RMFANOVA_DEFAULT, c("Cn_P1", "Dn_P1", "En_P1"))
  # the two that rejected true nulls at ~3x nominal must never be in the default
  expect_length(intersect(FCK_RMFANOVA_DEFAULT, FCK_RMFANOVA_INFLATED), 0)
  expect_length(intersect(FCK_RMFANOVA_DEFAULT, FCK_RMFANOVA_NOPOWER), 0)
  expect_setequal(FCK_RMFANOVA_INFLATED, c("En_P2", "En_B2"))
  expect_setequal(FCK_RMFANOVA_NOPOWER, c("Cn_P2", "Cn_B2"))
})

test_that("the global test is optional and reports why it did not run", {
  src <- code_of("server/50_fanova.R")
  expect_true(grepl("isTRUE(input$rm_global_test)", src, fixed = TRUE))
  expect_true(grepl("rm_global = rm_global", src, fixed = TRUE))
  ui <- paste(readLines(file.path(app_dir, "ui/50_fanova.R"), warn = FALSE), collapse = "\n")
  expect_true(grepl("complete balanced design", ui, fixed = TRUE))
  expect_true(grepl("anti-conservative", ui, fixed = TRUE))
  # absent package must be a reason, not an error
  source(file.path(app_dir, "server/09b_helpers_rmfanova.R"), local = TRUE)
  r <- fck_rmfanova_global(matrix(1, 5, 2), factor(c(1, 2)), factor(c("a", "a")))
  expect_false(r$ok)
  expect_true(nzchar(r$reason))
})

test_that("the deparsed kernels re-parse into valid functions", {
  # deparse() is only a safe way to export an estimator if the text it produces
  # is still that estimator. A long function can hit deparse's line cutoff and
  # emit something that no longer parses, which would break the exported script
  # silently -- the app would still be right, the download would not.
  src <- readLines(file.path(app_dir, "tests/real_data_run.R"), warn = FALSE)
  st <- grep("^extract_fns <- function", src)[1]
  en <- st - 1 + which(src[st:length(src)] == "}")[1]
  eval(parse(text = paste(src[st:en], collapse = "\n")))
  env <- extract_fns(file.path(app_dir, "server/50_fanova.R"),
                     c("perform_functional_anova", "perform_rm_fanova"))
  for (nm in c("perform_functional_anova", "perform_rm_fanova")) {
    f <- get(nm, envir = env)
    txt <- paste(deparse(f), collapse = "\n")
    g <- eval(parse(text = txt))
    expect_true(is.function(g), info = nm)
    expect_equal(names(formals(f)), names(formals(g)), info = nm)
    expect_gt(nchar(txt), 1000)      # it really is the whole function
  }
})
