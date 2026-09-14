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

# ==============================================================================
# The install offer (P21 phase 4)
# ==============================================================================
# An earlier audit removed automatic installation for three reasons, all sound:
# it silently changes the analysis, it cannot work where it is most needed, and
# it hides the real problem. Those argue against installing SILENTLY, not
# against installing at all. The offer must therefore keep every one of the
# three properties, and these checks are what stop a later convenience from
# quietly dropping one.

test_that("nothing is ever installed without an answer", {
  expect_true(grepl("dance_offer_install <- function", app, fixed = TRUE))
  fn <- regmatches(app, regexpr("dance_offer_install <- function(?s).*?\\n\\}", app, perl = TRUE))
  expect_length(fn, 1)
  # it must ASK, and must only install after a yes
  expect_true(grepl("readline(", fn, fixed = TRUE))
  expect_true(grepl('ans %in% c("y", "yes")', fn, fixed = TRUE))
  # the install call must come after that test, never before it
  expect_gt(regexpr("install.packages", fn, fixed = TRUE),
            regexpr('ans %in% c("y", "yes")', fn, fixed = TRUE))
})

test_that("it never even asks where an answer cannot be given", {
  fn <- regmatches(app, regexpr("dance_offer_install <- function(?s).*?\\n\\}", app, perl = TRUE))
  # reason 2 of the original audit: a deployed app, Rscript, CI or a locked-down
  # machine must fail immediately with the command, not pause trying to install
  expect_true(grepl("interactive()", fn, fixed = TRUE))
  expect_true(grepl("DANCE_NO_INSTALL", fn, fixed = TRUE))
  expect_true(grepl("if (!can_ask) return(pkgs)", fn, fixed = TRUE))
})

test_that("the packages are named before anything happens, and the versions after", {
  fn <- regmatches(app, regexpr("dance_offer_install <- function(?s).*?\\n\\}", app, perl = TRUE))
  # reason 3: name the missing package before acting on it
  expect_lt(regexpr('for \\(p in pkgs\\)', fn), regexpr("readline\\(", fn))
  # reason 1: what changed goes on the record
  expect_true(grepl("packageVersion", fn, fixed = TRUE))
  expect_true(grepl("Installed: ", fn, fixed = TRUE))
})

test_that("the required check still stops when the offer did not resolve it", {
  # the offer narrows the list; the stop must act on what is left, not on the
  # original -- otherwise a declined install starts the app anyway
  expect_true(grepl("missing_required <- dance_offer_install(missing_required)", app, fixed = TRUE))
  i_offer <- regexpr("missing_required <- dance_offer_install", app, fixed = TRUE)
  i_stop  <- regexpr("DANCE cannot start", app, fixed = TRUE)
  expect_gt(i_stop, i_offer)
})

test_that("the one-command installer reads its list from app.R", {
  inst <- paste(readLines(file.path(app_dir, "tools/install_dependencies.R"), warn = FALSE),
                collapse = "\n")
  # it must not repeat the package names -- a second copy is a second thing to
  # forget, which is the omission that started this
  expect_true(grepl("parse(file.path(root, \"app.R\"))", inst, fixed = TRUE))
  expect_false(grepl('"emmeans"', inst, fixed = TRUE))
  expect_false(grepl('"shinyWidgets"', inst, fixed = TRUE))
  # a missing REQUIRED package is a non-zero exit; a missing optional one is not
  expect_true(grepl("quit(status = 1)", inst, fixed = TRUE))
})
