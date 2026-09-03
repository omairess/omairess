# ==============================================================================
# tests/audit_test.R — run the audit regression suite WITHOUT testthat
#
# WHY THIS EXISTS
# ---------------
# tests/testthat/test-cosinor-audit.R is the deliverable the brief asked for.
# testthat could not be installed in the environment this work was done in
# (CRAN is blocked by the egress policy there), so the assertions in that file
# had never been executed -- and an unexecuted test is a claim, not a check.
#
# This runs the SAME file by supplying a minimal shim for the handful of
# testthat verbs it uses. It deliberately does not reimplement the assertions:
# duplicating them would let the two copies drift, and then the one that runs
# would stop being evidence about the one that ships.
#
# Run:  Rscript tests/audit_test.R          (from the FCK directory)
# ==============================================================================

app_dir <- if (dir.exists("server")) "." else "FCK"

.n_pass <- 0L; .n_fail <- 0L; .fails <- character(0); .ctx <- ""

.fail <- function(msg) {
  .n_fail <<- .n_fail + 1L
  .fails <<- c(.fails, sprintf("%s: %s", .ctx, msg))
}
.pass <- function() .n_pass <<- .n_pass + 1L

test_that <- function(desc, code) {
  .ctx <<- desc
  before_f <- .n_fail
  res <- tryCatch({ force(code); NULL },
                  error = function(e) conditionMessage(e))
  if (!is.null(res)) .fail(paste("ERROR:", res))
  if (.n_fail == before_f) cat(sprintf("ok   : %s\n", desc))
  else cat(sprintf("FAIL : %s\n", desc))
  invisible(NULL)
}

expect_equal <- function(object, expected, tolerance = 1e-8, ...) {
  cmp <- tryCatch(all.equal(object, expected, tolerance = tolerance,
                            check.attributes = FALSE),
                  error = function(e) conditionMessage(e))
  if (isTRUE(cmp)) .pass()
  else .fail(sprintf("expect_equal: %s", paste(cmp, collapse = "; ")))
  invisible(object)
}
expect_identical <- function(object, expected, ...) {
  if (identical(object, expected)) .pass() else .fail("expect_identical")
  invisible(object)
}
expect_true <- function(object, ...) {
  if (isTRUE(object)) .pass() else .fail("expect_true was not TRUE")
  invisible(object)
}
expect_false <- function(object, ...) {
  if (isFALSE(object)) .pass() else .fail("expect_false was not FALSE")
  invisible(object)
}
expect_null <- function(object, ...) {
  if (is.null(object)) .pass() else .fail("expect_null was not NULL")
  invisible(object)
}
expect_lt <- function(object, expected, ...) {
  if (isTRUE(object < expected)) .pass()
  else .fail(sprintf("expect_lt: %s !< %s", format(object), format(expected)))
  invisible(object)
}
expect_gt <- function(object, expected, ...) {
  if (isTRUE(object > expected)) .pass()
  else .fail(sprintf("expect_gt: %s !> %s", format(object), format(expected)))
  invisible(object)
}
expect_lte <- function(object, expected, ...) {
  if (isTRUE(object <= expected)) .pass()
  else .fail(sprintf("expect_lte: %s !<= %s", format(object), format(expected)))
  invisible(object)
}

cat("Running tests/testthat/test-cosinor-audit.R under a base-R shim\n")
cat(strrep("-", 68), "\n")
eval(parse(file.path(app_dir, "tests/testthat/test-cosinor-audit.R"),
           encoding = "UTF-8"),
     envir = environment())

cat(strrep("-", 68), "\n")
cat(sprintf("%d assertion(s) passed, %d failed.\n", .n_pass, .n_fail))
if (.n_fail > 0) {
  cat("\nFailures:\n"); for (f in .fails) cat("  -", f, "\n")
  quit(status = 1)
}
cat("Audit regression suite passed.\n")
