# ==============================================================================
# tests/testthat/test-fpca-varimax.R -- functional PCA with a VARIMAX rotation
#
# dance_fpca() is fda::pca.fd() and, on request, fda::varmx.pca.fd(); nothing
# is reimplemented. Checked here: that it IS those functions; that the rotation
# is a rotation (orthogonal, the curves the components rebuild unchanged, the
# total variance unchanged); and that fda's VARIMAX agrees with an independent
# implementation, stats::varimax (unnormalised), on the same component matrix.
# ==============================================================================
.libPaths(c("~/Rlib", .libPaths()))
`%||%` <- function(a, b) if (is.null(a)) b else a
app_dir <- if (dir.exists("server")) "." else if (dir.exists("../../server")) "../.." else "DANCE"
suppressMessages(library(fda))
ev <- new.env()
for (f in c("01_helpers_time.R", "03_helpers_clock.R", "04_helpers_fd.R"))
  sys.source(file.path(app_dir, "server", f), envir = ev)

mk_fd <- function(n = 30, seed = 3) {
  set.seed(seed); tt <- seq(0, 1, length.out = 25)
  Y <- t(sapply(seq_len(n), function(i) 2 * rnorm(1) * sin(2 * pi * tt) + rnorm(1) * cos(2 * pi * tt) +
                  1.5 * rnorm(1) * exp(-((tt - 0.7) / 0.08)^2) + rnorm(25, 0, 0.1)))
  smooth.basis(tt, t(Y), fdPar(create.bspline.basis(c(0, 1), 12), 2, 1e-6))$fd
}

test_that("no rotation is pca.fd itself", {
  fdo <- mk_fd()
  a <- ev$dance_fpca(fdo, 3, "none"); b <- pca.fd(fdo, nharm = 3)
  expect_identical(a$rotation, "none")
  expect_equal(a$scores, b$scores); expect_equal(a$varprop, b$varprop)
  expect_equal(a$harmonics$coefs, b$harmonics$coefs)
})

test_that("VARIMAX is fda::varmx.pca.fd, and it is a rotation", {
  fdo <- mk_fd()
  r <- ev$dance_fpca(fdo, 3, "varimax"); p0 <- pca.fd(fdo, nharm = 3); ref <- varmx.pca.fd(p0)
  expect_identical(r$rotation, "varimax")
  expect_equal(r$scores, ref$scores); expect_equal(r$harmonics$coefs, ref$harmonics$coefs)
  expect_equal(r$varprop_unrotated, p0$varprop)
  # orthogonal: R'R = I. To about 1e-7 only: fda's varmx() builds the rotation
  # from plane rotations with cos(45 deg) written as 0.7071068 (seven digits).
  expect_equal(crossprod(r$rotmat), diag(3), tolerance = 1e-6)
  # the same total variance, shared differently
  expect_equal(sum(r$varprop), sum(p0$varprop), tolerance = 1e-10)   # varmx.pca.fd rescales to it exactly
  expect_false(isTRUE(all.equal(r$varprop, p0$varprop)))
  # the curves the components rebuild are unchanged by the rotation
  x <- seq(0, 1, length.out = 101)
  rebuild <- function(p) eval.fd(x, p$harmonics) %*% t(p$scores)
  expect_equal(rebuild(r), rebuild(p0), tolerance = 1e-6)
})

test_that("fda's VARIMAX agrees with stats::varimax on the same component matrix", {
  fdo <- mk_fd(seed = 9)
  p0 <- pca.fd(fdo, nharm = 3)
  x <- seq(0, 1, length.out = 501)
  H <- eval.fd(x, p0$harmonics)
  R_fda <- varmx(H)
  R_st <- unclass(stats::varimax(H, normalize = FALSE, eps = 1e-12)$rotmat)
  # the same rotation up to the order and sign of its columns
  M <- abs(crossprod(R_fda, R_st))
  expect_equal(sort(apply(M, 1, max)), rep(1, 3), tolerance = 1e-3)
  # and the same criterion value (the variance of the squared loadings)
  crit <- function(L) sum(apply(L^2, 2, var))
  expect_equal(crit(H %*% R_fda), crit(H %*% R_st), tolerance = 1e-6)
})

test_that("one component is not rotated", {
  r <- ev$dance_fpca(mk_fd(), 1, "varimax")
  expect_identical(r$rotation, "none")
})
