# ==============================================================================
# tests/cosinor_convention_test.R
#
# ONE CONVENTION, ASSERTED END TO END.
#
# DANCE has had an acrophase units/convention discrepancy before, and the reason
# that class of bug survives review is that every individual line looks right:
# atan2 takes two arguments and both orders read plausibly, a sign flip is one
# character, and radians-to-hours differs from hours-to-radians by a division
# that is easy to write twice. Reading the code cannot settle it. Planting a
# KNOWN phase and demanding that every path return that same number can.
#
# The convention this file pins, everywhere in the app:
#
#   the model is        y = M + a cos(2 pi k t / P) + b sin(2 pi k t / P)
#   the acrophase is    phi = atan2(b, a), taken modulo 2 pi, so phi is in [0, 2pi)
#   phi is the time of the PEAK, positive, not the negative angle Bingham writes
#   in clock units       phi * P / (2 pi) / k
#   the kth harmonic has effective period P / k, so its phase wraps on P / k
#
# Each check below drives a different entry point with the same planted signal.
# If any one of them disagrees, one of them is wrong, and the failure names it.
# ==============================================================================

`%||%` <- function(a, b) if (is.null(a)) b else a
app <- if (dir.exists("server")) "." else ".."
e <- new.env()
for (f in c("07_helpers_circular", "08_helpers_cosinor", "08b_helpers_popcosinor",
            "08c_helpers_circstat", "08d_helpers_traj", "08e_helpers_trajfit",
            "08f_helpers_trajinf", "08g_helpers_trajband"))
  sys.source(file.path(app, "server", paste0(f, ".R")), envir = e)

pass <- 0L; fail <- 0L
chk <- function(ok, a, b) { if (isTRUE(ok)) { pass <<- pass + 1L; cat("ok   ", a, "\n") }
                            else { fail <<- fail + 1L; cat("FAIL ", b, "\n") } }

P <- 24

cat("-- the definition itself ---------------------------------------------\n")
# A pure cosine peaking at t = 0 has phi = 0; one peaking at t = 6 has phi = pi/2.
for (peak in c(0, 3, 6, 11, 18, 23)) {
  a <- cos(2 * pi * peak / P); b <- sin(2 * pi * peak / P)
  phi <- atan2(b, a) %% (2 * pi)
  hr <- e$phi_to_hours(phi, P, 1)
  chk(abs(((hr - peak + P/2) %% P) - P/2) < 1e-9,
      sprintf("a cosine peaking at t = %4.1f h gives phi = %.4f rad = %5.2f h", peak, phi, hr),
      sprintf("a cosine peaking at t = %.1f h reported %.2f h", peak, hr))
}

cat("\n-- the round trip, and the harmonic scaling --------------------------\n")
chk(max(abs(vapply(seq(0, 23.5, by = .5), function(h)
      e$phi_to_hours(e$hours_to_phi(h, P, 1), P, 1) - h, numeric(1)))) < 1e-9,
    "hours -> radians -> hours is the identity on the first harmonic",
    "the radians/hours round trip does not return what it was given")
# the kth harmonic lives on P/k: 2 h on harmonic 2 is a QUARTER of the way round
chk(abs(e$hours_to_phi(3, P, 2) - pi/2) < 1e-12 &&
      abs(e$phi_to_hours(pi/2, P, 2) - 3) < 1e-12,
    "the 2nd harmonic wraps on P/2 = 12 h, so 3 h is a quarter turn",
    "the harmonic scaling is not P/k")

cat("\n-- the two-stage path: fit_cosinor on a planted signal ---------------\n")
# fit_cosinor still lives inside the server body of 72_harmonic.R and cannot be
# sourced. It is SLICED OUT and run here, the same way phase 0 verified the
# circular helpers: comparing the shipped convention against a reimplementation
# would only test the reimplementation.
h72 <- readLines(file.path(app, "server", "72_harmonic.R"), warn = FALSE)
i0 <- grep("^  fit_cosinor <- function", h72)[1]
stopifnot(!is.na(i0))
depth <- 0L; i1 <- i0; started <- FALSE
repeat {
  depth <- depth + lengths(regmatches(h72[i1], gregexpr("\\{", h72[i1]))) -
                   lengths(regmatches(h72[i1], gregexpr("\\}", h72[i1])))
  if (depth > 0L) started <- TRUE
  if (started && depth <= 0L) break
  i1 <- i1 + 1L
  if (i1 > length(h72)) stop("fit_cosinor's closing brace was not found")
}
eval(parse(text = paste(h72[i0:i1], collapse = "\n")), envir = e)
cat(sprintf("     (fit_cosinor sliced from 72_harmonic.R, lines %d-%d)\n", i0, i1))

tp <- seq(0, 23, by = 1)
peak <- 8; amp <- 5
y <- 10 + amp * cos(2 * pi * (tp - peak) / P)
f1 <- suppressWarnings(e$fit_cosinor(tp, y, period = P, n_harmonics = 1,
                                     trend_type = "none"))
# The REPORTED fields, not the raw coefficients: what the app prints is what has
# to be right, and a convention error between the fit and the table would be
# invisible to a check on the coefficients alone.
phi_fit <- f1$acrophases[1] %% (2 * pi)
co <- f1$coefs   # c(mesor, trend..., cos1, sin1, ...)
chk(abs(e$phi_to_hours(phi_fit, P, 1) - peak) < 1e-6,
    sprintf("fit_cosinor reports the planted 8 h peak as %.4f h", e$phi_to_hours(phi_fit, P, 1)),
    sprintf("fit_cosinor put the peak at %.4f h, not 8", e$phi_to_hours(phi_fit, P, 1)))
chk(abs(f1$acrophases_time[1] - peak) < 1e-6,
    sprintf("and its own radians->hours column agrees (%.4f h)", f1$acrophases_time[1]),
    sprintf("acrophases_time says %.4f h while acrophases says %.4f h -- two conventions in one object",
            f1$acrophases_time[1], e$phi_to_hours(phi_fit, P, 1)))
chk(abs(f1$amplitudes[1] - amp) < 1e-6,
    "and the amplitude is the planted 5",
    "fit_cosinor's amplitude does not match the planted one")
chk(abs(atan2(co[3], co[2]) %% (2 * pi) - phi_fit) < 1e-9,
    "the reported acrophase is atan2(sin coefficient, cos coefficient), in that order",
    "the reported acrophase is not atan2(b, a) of its own coefficients")

cat("\n-- the one-stage path: the trajectory model, same planted signal -----\n")
n <- 10
Y <- t(vapply(seq_len(n), function(i)
  10 + amp * cos(2 * pi * (tp - peak) / P) + rnorm(length(tp), 0, .01), numeric(length(tp))))
set.seed(4)
d <- e$dance_traj_long(Y, tp, sprintf("S%02d", seq_len(n)),
                       list(Group = rep(c("a", "b"), each = n / 2)))
fit <- suppressWarnings(suppressMessages(
  e$dance_traj_fit(e$dance_traj_spec(d, P, 1, "none"))))
cf <- e$dance_traj_cell_coefs(fit, 1)
phi_traj <- atan2(cf$b[1], cf$a[1]) %% (2 * pi)
chk(abs(e$phi_to_hours(phi_traj, P, 1) - peak) < 1e-2,
    sprintf("the trajectory model recovers the same 8 h peak (%.4f h)",
            e$phi_to_hours(phi_traj, P, 1)),
    sprintf("the trajectory model put the peak at %.4f h", e$phi_to_hours(phi_traj, P, 1)))

# THE CROSS-CHECK THAT MATTERS: the two paths must agree to the noise, not
# merely each look reasonable on its own
chk(abs(e$phi_to_hours(phi_traj, P, 1) - e$phi_to_hours(phi_fit, P, 1)) < 1e-2,
    "and the two-stage and one-stage paths agree on the acrophase",
    "the two paths disagree on the acrophase -- one of them has the wrong convention")

cat("\n-- the derived-parameter tables use that same convention -------------\n")
ap_d <- e$dance_traj_amp_phase(cf)
ap_j <- e$dance_traj_amp_phase_joint(cf, n_draw = 4000)
chk(abs(ap_d$table$acrophase_rad[1] - phi_traj) < 1e-12 &&
      abs(ap_d$table$acrophase_time[1] - e$phi_to_hours(phi_traj, P, 1)) < 1e-9,
    "the delta-method table reports atan2(b, a) mod 2pi, converted with phi_to_hours",
    "the delta-method table's acrophase is not the one the coefficients imply")
chk(abs(ap_j$table$acrophase_rad[1] - ap_d$table$acrophase_rad[1]) < 1e-12,
    "and the joint-draw table reports the identical point estimate",
    "the two inference methods disagree on the POINT estimate, which cannot happen")
chk(identical(ap_d$effective_period, ap_j$effective_period) &&
      abs(ap_d$effective_period - P) < 1e-12,
    "both report the effective period P/k = 24 for the first harmonic",
    "the effective period differs between the two methods")

cat("\n-- the second harmonic wraps on 12 h, everywhere ---------------------\n")
y2 <- 10 + 3 * cos(2 * pi * 2 * (tp - 2) / P)
Y2 <- t(vapply(seq_len(n), function(i) y2 + rnorm(length(tp), 0, .01), numeric(length(tp))))
d2 <- e$dance_traj_long(Y2, tp, sprintf("S%02d", seq_len(n)),
                        list(Group = rep(c("a", "b"), each = n / 2)))
fit2 <- suppressWarnings(suppressMessages(
  e$dance_traj_fit(e$dance_traj_spec(d2, P, 2, "none"))))
cf2 <- e$dance_traj_cell_coefs(fit2, 2)
chk(abs(cf2$effective_period - 12) < 1e-12,
    "harmonic 2 declares an effective period of 12 h",
    sprintf("harmonic 2 reports an effective period of %.3f", cf2$effective_period))
ph2 <- atan2(cf2$b[1], cf2$a[1]) %% (2 * pi)
chk(abs(((ph2 * 12 / (2 * pi)) - 2 + 6) %% 12 - 6) < 5e-2,
    sprintf("and recovers the planted 2 h peak ON THAT 12 h period (%.3f h)",
            ph2 * 12 / (2 * pi)),
    sprintf("harmonic 2's peak came back at %.3f h on a 12 h period, not 2",
            ph2 * 12 / (2 * pi)))

cat("\n-- wrapping: a contrast is signed and never exceeds half the period --\n")
# two cells whose phases straddle the period boundary: the difference must be
# the SHORT way round, not the arithmetic difference of two clock readings
cf3 <- cf
cf3$a <- c(cos(2 * pi * 23 / P), cos(2 * pi * 1 / P))
cf3$b <- c(sin(2 * pi * 23 / P), sin(2 * pi * 1 / P))
cf3$cells <- c("late", "early"); cf3$joint <- list(diag(2) * 1e-6, diag(2) * 1e-6)
cf3$linfct <- NULL; cf3$vcov_beta <- NULL
pc <- e$dance_traj_phase_contrast(cf3, 1, 2)
chk(abs(pc$diff_time - 2) < 1e-6,
    sprintf("23 h vs 1 h is a contrast of +2 h, not -22 h (got %.4f)", pc$diff_time),
    sprintf("the phase contrast came back as %.4f h -- it took the long way round", pc$diff_time))
chk(abs(e$dance_traj_phase_contrast(cf3, 2, 1)$diff_time + 2) < 1e-6,
    "and reversing the pair reverses the sign",
    "the phase contrast is not antisymmetric in its arguments")

cat("\n-- the circular helpers share the convention -------------------------\n")
ang <- c(2 * pi * 23 / P, 2 * pi * 1 / P) * 1        # 23 h and 1 h as angles
m <- e$dance_circular_mean(ang) %% (2 * pi)
chk(abs(e$phi_to_hours(m, P, 1) %% P) < 1e-6 || abs(e$phi_to_hours(m, P, 1) - P) < 1e-6,
    sprintf("the circular mean of 23 h and 1 h is midnight (%.4f h), not noon",
            e$phi_to_hours(m, P, 1)),
    sprintf("the circular mean of 23 h and 1 h came back as %.4f h",
            e$phi_to_hours(m, P, 1)))

cat(sprintf("\n%s  (%d passed, %d failed)\n",
            if (fail == 0) "Cosinor convention tests PASSED" else "FAILURES", pass, fail))
if (fail > 0) quit(status = 1)
