# ==============================================================================
# tests/testthat/test-declared-packages.R
#
# Every package a server file reaches for must be DECLARED in app.R -- required
# if the app cannot run without it, optional if a feature cannot.
#
# The trajectory framework was built against emmeans, pbkrtest, lmerTest and
# glmmTMB without declaring any of them. The app therefore started cleanly on a
# machine that had none, drew the fitted curves, and gave the first sign of a
# problem as a terse refusal inside a panel the user had already waited two
# minutes for. Declaring them is what makes the startup check report them.
# ==============================================================================

app_dir <- if (dir.exists("server")) "." else if (dir.exists("../../server")) "../.." else "DANCE"
app <- paste(readLines(file.path(app_dir, "app.R"), warn = FALSE), collapse = "\n")
# COMMENT LINES ARE NOT CODE. A draft of this check scanned the raw text and
# reported `refund` as an undeclared dependency -- it appears only in a comment
# explaining that rmfanova declares refund but never calls it. The same mistake
# as the terminology check, which failed on the comment documenting the wording
# it was written to enforce.
srv_lines <- unlist(lapply(list.files(file.path(app_dir, "server"), pattern = "[.]R$",
                                      full.names = TRUE), readLines, warn = FALSE))
srv <- paste(grep("^\\s*#", srv_lines, value = TRUE, invert = TRUE), collapse = "\n")

# PARSE app.R, do not scrape it. A first draft pulled the two package vectors
# out with regular expressions, matched neither, and reported all sixteen
# declared packages as undeclared -- a check that fails on everything is as
# useless as one that fails on nothing. Evaluating the two assignments is exact.
declared_packages <- local({
  env <- new.env()
  for (ex in parse(file.path(app_dir, "app.R"))) {
    if (is.call(ex) && identical(as.character(ex[[1]]), "<-") &&
        is.name(ex[[2]]) &&
        as.character(ex[[2]]) %in% c("required_packages", "optional_packages"))
      eval(ex, env)
  }
  unique(c(env$required_packages, names(env$optional_packages)))
})

used_packages <- local({
  ns <- unlist(regmatches(srv, gregexpr('requireNamespace\\(\\s*"[A-Za-z0-9.]+"', srv)))
  ns <- gsub('.*"([A-Za-z0-9.]+)".*', "\\1", ns)
  qual <- unlist(regmatches(srv, gregexpr("\\b[A-Za-z][A-Za-z0-9.]*(?=::)", srv, perl = TRUE)))
  base_pkgs <- rownames(utils::installed.packages(priority = c("base", "recommended")))
  setdiff(unique(c(ns, qual)), c(base_pkgs, "utils", "stats", "grDevices", "graphics",
                                 "methods", "tools", "parallel", "compiler", "datasets"))
})

test_that("the two package vectors are found at all", {
  # the guard the first draft needed: if this is empty the check below passes
  # vacuously or fails on everything, and either way tells you nothing
  expect_gt(length(declared_packages), 10)
  expect_gt(length(used_packages), 3)
})

test_that("every package the server reaches for is declared in app.R", {
  missing <- setdiff(used_packages, declared_packages)
  expect_equal(missing, character(0),
               info = paste("undeclared:", paste(missing, collapse = ", ")))
})

test_that("the trajectory framework's four dependencies say what they cost", {
  opt <- regmatches(app, regexpr("optional_packages <- c\\((?s).*?\\n\\)", app, perl = TRUE))
  for (p in c("emmeans", "pbkrtest", "lmerTest", "glmmTMB"))
    expect_true(grepl(sprintf("\\n\\s*%s\\s*=\\s*\"[^\"]{10,}\"", p), opt), info = p)
  # and the consequence is named, not just the feature
  expect_true(grepl("falls back to Satterthwaite", opt, fixed = TRUE))
  expect_true(grepl("falls back to a likelihood-ratio test", opt, fixed = TRUE))
})

test_that("the emmeans refusal tells the user what to do and what still works", {
  inf <- paste(readLines(file.path(app_dir, "server/08f_helpers_trajinf.R"), warn = FALSE),
               collapse = "\n")
  expect_true(grepl('install.packages(\\"emmeans\\")', inf, fixed = TRUE))
  expect_true(grepl("are unaffected", inf, fixed = TRUE))
})
