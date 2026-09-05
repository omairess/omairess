# ==============================================================================
# tests/warping_test.R -- the registration estimator, executed
#
# The review found that time-warped PCA added runif() noise to an estimated
# transformation, attenuated the measured lag to 5% of itself, and in the
# landmark branch returned the curves unwarped while the UI announced success.
# Source greps cannot show that a rewrite fixed it; only running it can.
#
#   Rscript tests/warping_test.R
# ==============================================================================
.libPaths(c("~/Rlib", .libPaths()))
suppressWarnings(suppressMessages(library(fda)))
app <- if (dir.exists("server")) "." else "FCK"
# reuse the AST-based extractor the other harnesses use
src_txt <- readLines(file.path(app, "tests/real_data_run.R"), warn = FALSE)
st <- grep("^extract_fns <- function", src_txt)[1]
en <- st - 1 + which(src_txt[st:length(src_txt)] == "}")[1]
eval(parse(text = paste(src_txt[st:en], collapse = "\n")))
env <- extract_fns(file.path(app, "server/40_fpca.R"),
                   c("linear_shift_alignment", "parametric_alignment"))
env$linear_shift_alignment <- get("linear_shift_alignment", envir = env)
env$parametric_alignment  <- get("parametric_alignment",  envir = env)

# curves with a KNOWN phase shift
set.seed(9)
tp <- seq(0, 1, length.out = 60); n <- 24
true_shift <- runif(n, -0.15, 0.15)
Y <- sapply(1:n, function(i) sin(2*pi*(tp - true_shift[i])) + rnorm(60, 0, .05))
b  <- create.bspline.basis(c(0,1), nbasis = 20)
fd <- smooth.basis(tp, Y, fdPar(b, 2, 1e-6))$fd

cat("=== DETERMINISM ===\n")
r1 <- env$linear_shift_alignment(fd, time_points = tp)
r2 <- env$linear_shift_alignment(fd, time_points = tp)
cat(sprintf("  identical warp functions across two runs: %s\n",
            isTRUE(all.equal(r1$warp_functions, r2$warp_functions))))
cat(sprintf("  identical registered curves            : %s\n",
            isTRUE(all.equal(r1$registered_curves, r2$registered_curves))))

cat("\n=== MONOTONICITY (every warp must be strictly increasing) ===\n")
mono <- apply(r1$warp_functions, 2, function(w) all(diff(w) > 0))
cat(sprintf("  monotone warps: %d of %d\n", sum(mono), length(mono)))

cat("\n=== RECOVERY OF THE KNOWN SHIFT ===\n")
est <- r1$shifts
cat(sprintf("  correlation(estimated, true) = %.3f\n", cor(est, true_shift)))
cat(sprintf("  slope of est on true         = %.3f  (1.0 = correct scale)\n",
            coef(lm(est ~ true_shift))[2]))
cat(sprintf("  old code's slope would be    = %.3f  (0.1 x 0.5 attenuation)\n",
            coef(lm(est ~ true_shift))[2] * 0.05))

cat("\n=== IDENTITY IN -> IDENTITY OUT ===\n")
Yc <- sapply(1:8, function(i) sin(2*pi*tp))
fdc <- smooth.basis(tp, Yc, fdPar(b, 2, 1e-6))$fd
rc <- env$linear_shift_alignment(fdc, time_points = tp)
cat(sprintf("  max |shift| on identical curves: %.4g\n", max(abs(rc$shifts))))
cat(sprintf("  max |warp - t|                 : %.4g\n",
            max(abs(rc$warp_functions - tp))))

cat("\n=== PARAMETRIC: monotone over the restricted range ===\n")
rp <- env$parametric_alignment(fd, family = "quadratic", param_range = c(0.5, 2),
                               time_points = tp)
mp <- apply(rp$warp_functions, 2, function(w) all(diff(w) > -1e-12))
cat(sprintf("  monotone warps: %d of %d ; alpha range used [%.3f, %.3f]\n",
            sum(mp), length(mp), min(rp$alpha_values), max(rp$alpha_values)))

# ---- assertions -------------------------------------------------------------
stopifnot(
  "warping must be deterministic"     = isTRUE(all.equal(r1$warp_functions, r2$warp_functions)),
  "registration must be deterministic"= isTRUE(all.equal(r1$registered_curves, r2$registered_curves)),
  "every warp must be strictly increasing" = all(mono),
  "the estimated shift must track the true one" = cor(est, true_shift) > 0.9,
  "the shift must not be attenuated"  = coef(lm(est ~ true_shift))[2] > 0.7,
  "identical curves must give the identity warp" = max(abs(rc$warp_functions - tp)) < 1e-8,
  "the parametric family must stay monotone" = all(mp)
)
cat("\nWarping tests passed.\n")
