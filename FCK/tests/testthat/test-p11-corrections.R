# ==============================================================================
# tests/testthat/test-p11-corrections.R — the seventh review's findings
#
#   P11.1  server/90_export.R wrote its own copy of the registration algorithms
#          with add("..."), and the copy predated P0.8. The exported script
#          recovered about a twelfth of each periodic shift AND REVERSED ITS
#          SIGN; the parametric branch was the pre-P4.1 code and deformed the
#          time axis by 0.4999 of the domain on curves needing no warping; the
#          `logistic` family was missing from its switch; and the landmark
#          branch emitted a comment and no code, leaving `registered_curves`
#          undefined. So "Export R code" reproduced the analysis for none of the
#          three methods.
#
#   P11.2  found while fixing P11.1, and worse, because it affected the LIVE
#          app: fck_landmark_warp() was defined at depth 11 inside the landmark
#          branch of a reactive observer, while its only caller was defined at
#          the top level of the same file. R scopes functions lexically, so the
#          call could never resolve; the surrounding tryCatch turned the error
#          into NULL, and the caller treated NULL as "warping failed" and
#          substituted a LINEAR SHIFT. Every landmark registration silently
#          returned a different method's answer.
#
#   P11.3  the smoothing diagnostics still hand-built their basis over
#          rangeval = c(1, n_time) with the column index as argvals, while
#          production uses elapsed hours when real-time smoothing is on.
#
# The structural fix for all three is the same one this audit has now reached
# five times: ONE definition, in a pure file, shared by the app, the tests and
# the exported script. See tests/warp_export_roundtrip_test.R and
# tests/diagnostic_axis_test.R for the executable halves.
# ==============================================================================
.libPaths(c("~/Rlib", .libPaths()))
app_dir <- if (dir.exists("server")) "." else if (dir.exists("../../server")) "../.." else "FCK"

code_of <- function(f) {
  ln <- readLines(file.path(app_dir, f), warn = FALSE)
  ln <- sub("#.*$", "", ln)
  paste(ln, collapse = "\n")
}
kern_env <- function() {
  e <- new.env(parent = globalenv())
  source(file.path(app_dir, "server/05_helpers_warp.R"), local = e)
  e
}

# ======================================================== P11.2 the live bug ==
test_that("P11.2: every registration kernel is a top-level, reachable definition", {
  e <- kern_env()
  for (f in c("fck_landmark_warp", "linear_shift_alignment",
              "parametric_alignment", "landmark_alignment_simple"))
    expect_true(exists(f, envir = e, inherits = FALSE), info = f)
  # the defect precisely: the caller's enclosing environment must contain the
  # helper. While fck_landmark_warp lived in an observer's frame it did not.
  expect_identical(environment(e$landmark_alignment_simple),
                   environment(e$fck_landmark_warp))
})

test_that("P11.2: landmark registration returns a result instead of NULL", {
  skip_if_not_installed("fda")
  suppressMessages(library(fda))
  e  <- kern_env()
  tp <- seq(0, 1, length.out = 60)
  Y  <- sapply(1:5, function(i) sin(2*pi*(tp - i*0.01)) + 0.3*sin(4*pi*tp))
  fd <- smooth.basis(tp, Y, fdPar(create.bspline.basis(c(0,1), nbasis = 15), 2, 1e-8))$fd

  auto <- e$landmark_alignment_simple(fd, c(0.25, 0.5, 0.75), tp)
  expect_false(is.null(auto))
  expect_identical(auto$method, "landmark")          # NOT the shift fallback
  expect_true(all(apply(auto$warp_functions, 2, function(h) all(diff(h) > 0))))

  man <- e$landmark_alignment_simple(fd, c(0.25, 0.5, 0.75), tp,
                                     landmark_points = data.frame(x = c(0.3, 0.6)))
  expect_false(is.null(man))
  expect_identical(man$method, "landmark")
  expect_equal(man$landmarks_used, c(0.3, 0.6))
})

test_that("P11.2: the kernels are pure, so they can be deparsed into a script", {
  skip_if_not_installed("codetools")
  e <- kern_env()
  for (f in c("fck_landmark_warp", "linear_shift_alignment",
              "parametric_alignment", "landmark_alignment_simple")) {
    g <- codetools::findGlobals(get(f, envir = e), merge = TRUE)
    expect_identical(intersect(g, c("values", "input", "session",
                                    "showNotification")), character(0), info = f)
  }
  # the rejection notice is returned, not raised, so the caller decides
  expect_true(grepl("rejection_warning = rejection_warning",
                    code_of("server/05_helpers_warp.R"), fixed = TRUE))
  expect_true(grepl("showNotification(warp_fd$rejection_warning",
                    code_of("server/40_fpca.R"), fixed = TRUE))
})

# ================================================== P11.1 the stale export ====
test_that("P11.1: the generator emits the kernels instead of rewriting them", {
  ex <- code_of("server/90_export.R")
  for (k in c("linear_shift_alignment", "parametric_alignment",
              "landmark_alignment_simple"))
    expect_true(grepl(paste0('emit_kernel("', k, '")'), ex, fixed = TRUE), info = k)
  # the signatures of the pre-P0.8 copy, matched as code the generator would
  # WRITE. code_of() strips comments, so the AUDIT note explaining the removal
  # cannot satisfy the guard -- a mistake made three times in this audit.
  expect_false(grepl('shifts[i] <- best_lag / n_time * 0.1', ex, fixed = TRUE))
  expect_false(grepl('time_points - shifts[i] * 0.5', ex, fixed = TRUE))
  expect_false(grepl('warp_functions[1,i] <- 0', ex, fixed = TRUE))
  expect_false(grepl("if(abs(alpha-1) < 0.001) t", ex, fixed = TRUE))
  # the landmark branch used to emit nothing at all, leaving the variable the
  # next section reads undefined
  expect_false(grepl("(Implementation depends on detected landmarks)", ex, fixed = TRUE))
})

test_that("P11.1: the export no longer asks for statistics that were removed", {
  ex <- code_of("server/90_export.R")
  # AIC/BIC went at P5.2 and "variance explained by warping" at P5.3. add()
  # pastes its arguments, so sprintf() of a NULL field is character(0) and the
  # whole line silently vanished rather than erroring.
  expect_false(grepl("stats$AIC", ex, fixed = TRUE))
  expect_false(grepl("stats$BIC", ex, fixed = TRUE))
  expect_false(grepl("variance_explained_by_warping", ex, fixed = TRUE))
  # fields are read exactly, not by $ partial matching (the P6.2 lesson)
  expect_true(grepl('v <- stats[[key]]', ex, fixed = TRUE))
  # and the second basis rule is gone: the kernel already returns regfd
  expect_false(grepl("reg_basis <- create.bspline.basis(c(0, 1), nbasis = min(20, n_time - 2))",
                     ex, fixed = TRUE))
  expect_true(grepl('add("reg_fd <- warp_result$regfd")', ex, fixed = TRUE))
})

test_that("P11.1: the exported kernel reproduces the live one exactly", {
  skip_if_not_installed("fda")
  skip_if_not_installed("codetools")
  suppressMessages(library(fda))
  e  <- kern_env()
  tp <- seq(0, 1, length.out = 100)
  prof <- function(t) sin(2*pi*t) + 0.4*sin(4*pi*t + 0.8)
  Y  <- sapply(c(0, 0.05, 0.10, -0.05, -0.10), function(d) prof(tp - d))
  fd <- smooth.basis(tp, Y, fdPar(create.fourier.basis(c(0,1), nbasis = 15), 2, 1e-8))$fd

  # deparse() then re-evaluate, which is exactly what emit_kernel() does
  txt <- paste(deparse(e$linear_shift_alignment), collapse = "\n")
  ev  <- new.env(parent = globalenv())
  eval(parse(text = paste("linear_shift_alignment <-", txt)), envir = ev)

  live <- e$linear_shift_alignment(fd, periodic = TRUE, reference = "first", time_points = tp)
  expd <- ev$linear_shift_alignment(fd, periodic = TRUE, reference = "first", time_points = tp)
  expect_equal(live$shifts, expd$shifts)
  expect_equal(live$registered_curves, expd$registered_curves)
  # and the magnitude is the ESTIMATE, not a twelfth of it
  expect_equal(abs(live$shifts[3]), 0.10, tolerance = 0.02)
})

test_that("P11.1: every family the UI offers is reachable and leaves identity alone", {
  skip_if_not_installed("fda")
  suppressMessages(library(fda))
  e  <- kern_env()
  tp <- seq(0, 1, length.out = 100)
  prof <- function(t) sin(2*pi*t) + 0.4*sin(4*pi*t + 0.8)
  Z  <- sapply(1:5, function(i) prof(tp))
  fz <- smooth.basis(tp, Z, fdPar(create.fourier.basis(c(0,1), nbasis = 15), 2, 1e-8))$fd
  ui <- code_of("ui/40_settings.R")
  for (fam in c("power", "exponential", "quadratic", "logistic")) {
    expect_true(grepl(paste0('"', fam, '"'), ui, fixed = TRUE), info = fam)
    r <- e$parametric_alignment(fz, family = fam, param_range = c(0.5, 2), time_points = tp)
    expect_lt(max(abs(r$warp_functions - tp)), 0.02)
    expect_gte(r$identity_alpha, r$param_range_used[1])
    expect_lte(r$identity_alpha, r$param_range_used[2])
  }
})

# ============================================= P11.3 the diagnostic axis ======
test_that("P11.3: the diagnostics build no axis or production basis of their own", {
  for (f in c("server/20_smoothing.R", "server/30_diagnostics.R")) {
    src <- code_of(f)
    expect_false(grepl("create.fourier.basis(rangeval = c(1, n_time)", src, fixed = TRUE), info = f)
    expect_false(grepl("create.bspline.basis(rangeval = c(1, n_time)", src, fixed = TRUE), info = f)
    expect_false(grepl("using_real_time <- !is.null(real_time)", src, fixed = TRUE), info = f)
  }
  expect_true(grepl("axis            <- fck_smoothing_axis(input, values)",
                    code_of("server/20_smoothing.R"), fixed = TRUE))
  dg <- code_of("server/30_diagnostics.R")
  expect_true(grepl("axis  <- fck_smoothing_axis(input, values)", dg, fixed = TRUE))
  expect_true(grepl("axis   <- fck_smoothing_axis(input, values)", dg, fixed = TRUE))
  # argvals must be the axis the basis was built over, in both sweep branches
  expect_true(grepl("fck_auto_lambda(data_mat, axis$t_full, basis", dg, fixed = TRUE))
  expect_true(grepl("smooth.basis(axis$t_full[valid], y[valid], fdP)", dg, fixed = TRUE))
})

test_that("P11.3: production and diagnostic bases agree on every configuration", {
  skip_if_not_installed("fda")
  suppressMessages(library(fda))
  e <- new.env(parent = globalenv())
  for (f in c("server/01_helpers_time.R", "server/03_helpers_clock.R",
              "server/04_helpers_fd.R")) source(file.path(app_dir, f), local = e)

  hrs    <- c(0,1,2,3,4,5,6,8,10,12,15,18,21,23)   # unevenly spaced on purpose
  labels <- sprintf("%02d:00", hrs)
  n_time <- length(labels)
  vals   <- list(data = matrix(0, 5, n_time), time_labels = labels)

  for (real in c(FALSE, TRUE)) for (cyc in c(FALSE, TRUE)) {
    inp  <- list(use_real_time = real, is_cyclic = cyc, smooth_method = "manual",
                 n_basis = 12, n_basis_manual = 12)
    axis <- e$fck_smoothing_axis(inp, vals)
    b    <- e$fck_smoothing_basis(axis, e$fck_smoothing_nbasis(inp, n_time), "manual")
    expect_equal(b$rangeval, range(axis$t_full))
  }

  # the divergence P11.3 measured: with real time on, production spans elapsed
  # hours [0, 23] and, when cyclic, uses period 24 -- where the hand-built
  # diagnostic used [1, 14] and period 13, a 13-hour rhythm fitted to a day.
  inp  <- list(use_real_time = TRUE, is_cyclic = TRUE, smooth_method = "manual",
               n_basis = 12, n_basis_manual = 12)
  axis <- e$fck_smoothing_axis(inp, vals)
  expect_true(axis$using_real_time)
  b <- e$fck_smoothing_basis(axis, e$fck_smoothing_nbasis(inp, n_time), "manual")
  expect_equal(b$rangeval, c(0, 23))
  expect_equal(as.numeric(b$params), 24)

  old <- create.fourier.basis(rangeval = c(1, n_time), nbasis = min(n_time, 13))
  expect_false(isTRUE(all.equal(old$rangeval, b$rangeval)))
})

test_that("P11.3: nb_fourier varies only the count, and only for the sweep", {
  skip_if_not_installed("fda")
  suppressMessages(library(fda))
  e <- new.env(parent = globalenv())
  for (f in c("server/01_helpers_time.R", "server/03_helpers_clock.R",
              "server/04_helpers_fd.R")) source(file.path(app_dir, f), local = e)
  vals <- list(data = matrix(0, 5, 14), time_labels = sprintf("%02d:00", c(0:6,8,10,12,15,18,21,23)))
  inp  <- list(use_real_time = TRUE, is_cyclic = TRUE, smooth_method = "manual",
               n_basis = 12, n_basis_manual = 12)
  axis <- e$fck_smoothing_axis(inp, vals)
  prod <- e$fck_smoothing_basis(axis, 12, "manual")
  expect_equal(e$fck_smoothing_basis(axis, 12, "manual", nb_fourier = NULL)$nbasis,
               prod$nbasis)
  sw <- e$fck_smoothing_basis(axis, 12, "manual", nb_fourier = 7)
  expect_equal(sw$nbasis, 7)
  expect_equal(sw$rangeval, prod$rangeval)
  expect_equal(sw$params, prod$params)
})

# ====================================== P11.4 the resolution of the estimate ==
test_that("P11.4: the phase estimate reports its own resolution", {
  src <- code_of("server/40_fpca.R")
  expect_true(grepl("Resolution of the phase estimate", src, fixed = TRUE))
  expect_true(grepl("1 / (n_grid - 1)", src, fixed = TRUE))
  # the grid must NOT have been silently densified: registration still runs on
  # the analysis grid, so no previously reported number moves.
  expect_true(grepl("time_grid <- seq(0, 1, length.out = n)", src, fixed = TRUE))
  expect_false(grepl("length.out = 256", src, fixed = TRUE))
})

test_that("P11.4: the reported resolution is the estimator's actual floor", {
  skip_if_not_installed("fda")
  suppressMessages(library(fda))
  e <- kern_env()
  # A displacement well below one grid step cannot be resolved; one of several
  # steps can. This is the property the on-screen note describes.
  for (n in c(24L, 49L)) {
    tp   <- seq(0, 1, length.out = n)
    step <- 1 / (n - 1)
    prof <- function(t) exp(cos(2 * pi * t))
    Y  <- sapply(c(0, 0.2 * step, 3 * step), function(d) prof(tp - d))
    fd <- smooth.basis(tp, Y, fdPar(create.fourier.basis(c(0,1), nbasis = min(n, 13)), 2, 1e-8))$fd
    r  <- e$linear_shift_alignment(fd, periodic = TRUE, reference = "first", time_points = tp)
    # every estimate lands on the lattice of grid steps
    expect_equal(r$shifts / step, round(r$shifts / step), tolerance = 1e-6,
                 info = paste("n =", n))
    # a sub-step displacement is reported as zero, not as a small real shift
    expect_equal(r$shifts[2], 0, info = paste("n =", n))
    # a displacement of several steps is recovered
    expect_equal(abs(r$shifts[3]), 3 * step, tolerance = step, info = paste("n =", n))
  }
})

test_that("P11.4: the publication report states the resolution of the shift", {
  src <- code_of("server/93_apa_report.R")
  expect_true(grepl("resolution of the estimate is one grid step", src, fixed = TRUE))
  expect_true(grepl("differences smaller than that are not", src, fixed = TRUE))
  # stated in the Methods, next to the method, not buried in the caveats
  expect_lt(regexpr("resolution of the estimate is one grid step", src, fixed = TRUE),
            regexpr("Registration diagnostics", src, fixed = TRUE))
})
