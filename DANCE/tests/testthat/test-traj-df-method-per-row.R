# ==============================================================================
# tests/testthat/test-traj-df-method-per-row.R
#
# A user asked why the denominator df of one block was 1300 and of the next
# 2443, on the same fit. Most of that gap is real and statistical -- a block
# whose contrasts are purely between-participant sits at the between-participant
# df, and a block mixing in directions the residual informs is pulled above it.
# But dance_traj_block_test() asks for Kenward-Roger and falls through to
# Satterthwaite when KRmodcomp fails, SILENTLY and one block at a time, and the
# table read the method off the FIRST successful result and printed it under
# every row. So a table could carry two approximations and claim one, and a df
# that is merely from the other approximation looks like a bug in the fit.
#
# These assert the PROPERTY on the parsed source, not its formatting: the file
# is re-deparsed first, so comments and line breaks cannot satisfy or break a
# check. (The renderers live inside the monolithic server function and cannot
# be called in isolation.)
# ==============================================================================
.libPaths(c("~/Rlib", .libPaths()))
app_dir <- if (dir.exists("server")) "." else if (dir.exists("../../server")) "../.." else "DANCE"

# canonical form: parsed and re-deparsed, all runs of whitespace collapsed
canon <- function(f) {
  x <- paste(vapply(parse(file.path(app_dir, f)),
                    function(e) paste(deparse(e), collapse = " "), character(1)),
             collapse = " ")
  gsub("[[:space:]]+", " ", x)
}

test_that("the table's method claim is taken over ALL tests, not the first one", {
  s <- canon("server/72_harmonic.R")
  # methods_all collects every distinct method actually used
  expect_true(grepl("methods_all = unique(vapply(Filter(function(r) isTRUE(r$ok), computed)",
                    s, fixed = TRUE))
  # and the per-row marker is driven by whether that set has more than one entry
  expect_true(grepl("mixed <- length(tt$methods_all %||% character(0)) > 1L",
                    s, fixed = TRUE))
})

test_that("each block row can show the method that produced it", {
  s <- canon("server/72_harmonic.R")
  # the row reads its OWN result's method ...
  expect_true(grepl("short_m(r$method", s, fixed = TRUE))
  # ... and short_m distinguishes the approximations that differ in ddf
  expect_true(grepl("Kenward", s, fixed = TRUE))
  expect_true(grepl("Satterth", s, fixed = TRUE))
})

test_that("short_m maps every method string dance_traj_block_test can return", {
  # Pull the helper out of the source and exercise it, so the mapping is tested
  # rather than merely present.
  src <- paste(readLines(file.path(app_dir, "server/72_harmonic.R"), warn = FALSE),
               collapse = "\n")
  i <- regexpr("short_m <- function(m)", src, fixed = TRUE)
  expect_gt(i, 0)
  tail_src <- substr(src, i, nchar(src))
  # the definition ends at the first line that starts a new statement at 4 spaces
  j <- regexpr("\n    rows <- ", tail_src, fixed = TRUE)
  expect_gt(j, 0)
  short_m <- eval(parse(text = substr(tail_src, 1, j)))
  expect_equal(short_m("Kenward-Roger F"), "KR")
  expect_equal(short_m("Satterthwaite F"), "Satt.")
  expect_equal(short_m("likelihood-ratio on ML refits"), "LRT")
  expect_equal(short_m("asymptotic likelihood-ratio (glmmTMB)"), "Wald")
})

test_that("the fall-through that makes a table non-uniform still exists", {
  # If KR ever stops falling back, these markers become dead code rather than
  # wrong -- but the fall-through is what the marker exists to disclose, so a
  # silent removal of it should be noticed here.
  s <- canon("server/08f_helpers_trajinf.R")
  expect_true(grepl("pbkrtest::KRmodcomp(mL, m0)", s, fixed = TRUE))
  # KR returns only on success; failure falls past it to the Satterthwaite path
  i_kr   <- regexpr("KRmodcomp(mL, m0)", s, fixed = TRUE)
  i_satt <- regexpr("lmerTest::contest(m, L, joint = TRUE)", s, fixed = TRUE)
  expect_gt(i_satt, i_kr)
  # and each path labels itself, which is what the row marker reads
  expect_true(grepl('method = "Kenward-Roger F"', s, fixed = TRUE))
  expect_true(grepl('method = "Satterthwaite F"', s, fixed = TRUE))
})
