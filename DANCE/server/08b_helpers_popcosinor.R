# ==============================================================================
# server/08b_helpers_popcosinor.R — the population-mean cosinor, m groups
# ==============================================================================
# What the app has done until now is a TWO-STAGE comparison: fit a cosinor per
# participant, then compare the resulting point estimates between groups with a
# t test (or Watson-Williams on the acrophases). That is a legitimate procedure
# and it is what the Cosinor: pairwise tests tab does. It has two costs:
#
#   it treats each participant's estimate as if it were measured without error,
#   so a participant whose own rhythm is barely identified counts as much as one
#   whose rhythm is precise; and
#
#   it tests amplitude and acrophase SEPARATELY on quantities that are two
#   coordinates of one bivariate object, so the pair of tests is not a test of
#   "the rhythm differs".
#
# Bingham, Arbogast, Cornelissen-Guillaume, Lee & Halberg (1982), Chronobiologia
# 9(4), 397-439, give the population-mean cosinor instead: the per-participant
# (cosine, sine) coefficient pairs are averaged as VECTORS, the within-group
# covariance of that pair is pooled across groups, and MESOR, amplitude and
# acrophase are each tested as an F ratio against the part of that covariance
# that matters for them.
#
# The geometry is what makes the formulas checkable rather than transcribed. Let
# u = (cos phi, sin phi) be the unit vector along the pooled mean acrophase.
#
#   AMPLITUDE is the length of the mean vector, so its sampling variance is the
#   variance of the coefficient pair ALONG u:
#       s2_along = s2_beta cos^2 + 2 s_betagamma cos sin + s2_gamma sin^2
#
#   ACROPHASE is the angle, so its variance is the variance PERPENDICULAR to u,
#   scaled by the squared amplitude (a fixed angular error is a larger
#   displacement on a longer vector):
#       s2_perp  = s2_beta sin^2 - 2 s_betagamma cos sin + s2_gamma cos^2
#
# Bingham writes these with the opposite sign on the cross term because his
# acrophase is NEGATIVE by convention. Using atan2(gamma, beta) throughout, the
# signs above are the consistent ones -- and calibration under the null is what
# settles it, not the transcription: tests/pop_cosinor_test.R simulates each null
# and checks the rejection rate against the nominal level.
#
# THE CAVEAT THAT TRAVELS WITH THE TEST. Bingham et al. note that an amplitude
# difference cannot be interpreted when the groups also differ in acrophase: the
# amplitudes are then being compared about different phases. This module returns
# the acrophase verdict alongside the amplitude one so the caller can say so,
# rather than leaving the reader to remember.
#
# Pure: coefficients in, list out.
# ==============================================================================

if (!exists("%||%", mode = "function")) `%||%` <- function(a, b) if (is.null(a)) b else a

# beta  cosine coefficients, one per participant
# gamma sine coefficients,   one per participant
# mesor per-participant MESOR
# group a factor with m >= 2 levels
dance_pop_cosinor <- function(beta, gamma, mesor, group, period = 24, harmonic = 1) {
  keep <- is.finite(beta) & is.finite(gamma) & is.finite(mesor) & !is.na(group)
  beta <- beta[keep]; gamma <- gamma[keep]; mesor <- mesor[keep]
  group <- droplevels(factor(group[keep]))
  gl <- levels(group); m <- length(gl); K <- length(beta)

  if (m < 2) return(list(ok = FALSE, message = "The population-mean cosinor needs at least 2 groups."))
  kj <- as.integer(table(group))
  if (any(kj < 2))
    return(list(ok = FALSE, message = paste0(
      "Every group needs at least 2 participants for a within-group covariance; ",
      "group sizes are ", paste(sprintf("%s = %d", gl, kj), collapse = ", "), ".")))
  if (K - m < 1) return(list(ok = FALSE, message = "Too few participants for the error degrees of freedom."))

  df1 <- m - 1L; df2 <- K - m

  # ---- group means, as VECTORS ----------------------------------------------
  bbar <- tapply(beta,  group, mean)[gl]
  gbar <- tapply(gamma, group, mean)[gl]
  Mbar <- tapply(mesor, group, mean)[gl]
  Aj   <- sqrt(bbar^2 + gbar^2)
  phij <- atan2(gbar, bbar) %% (2 * pi)

  # ---- pooled within-group (co)variance of the coefficient pair -------------
  cb <- beta  - bbar[as.character(group)]
  cg <- gamma - gbar[as.character(group)]
  cm <- mesor - Mbar[as.character(group)]
  s2b  <- sum(cb^2)  / df2
  s2g  <- sum(cg^2)  / df2
  sbg  <- sum(cb*cg) / df2
  s2M  <- sum(cm^2)  / df2

  # ---- the pooled mean direction --------------------------------------------
  # A doubled-angle mean, so two acrophases half a cycle apart do not average to
  # something between them: this is a direction, not a number.
  phi_t <- 0.5 * atan2(sum(kj * Aj^2 * sin(2 * phij)),
                       sum(kj * Aj^2 * cos(2 * phij)))
  cs <- cos(phi_t); sn <- sin(phi_t)
  s2_along <- s2b * cs^2 + 2 * sbg * cs * sn + s2g * sn^2      # amplitude direction
  s2_perp  <- s2b * sn^2 - 2 * sbg * cs * sn + s2g * cs^2      # acrophase direction

  Fp <- function(num, den) {
    if (!is.finite(den) || den <= 0) return(c(NA_real_, NA_real_))
    f <- num / den
    c(f, stats::pf(f, df1, df2, lower.tail = FALSE))
  }
  Abar <- sum(kj * Aj) / K
  Mgrand <- sum(kj * Mbar) / K

  fM <- Fp(sum(kj * (Mbar - Mgrand)^2) / df1, s2M)
  fA <- Fp(sum(kj * (Aj - Abar)^2) / df1, s2_along)
  fP <- Fp(sum(kj * Aj^2 * sin(phij - phi_t)^2) / df1, s2_perp)

  acro_differs <- is.finite(fP[2]) && fP[2] < 0.05

  list(
    ok = TRUE, n_groups = m, n_subjects = K, group_names = gl, group_sizes = kj,
    df1 = df1, df2 = df2, period = period, harmonic = harmonic,
    groups = data.frame(
      group = gl, n = kj, mesor = as.numeric(Mbar),
      amplitude = as.numeric(Aj),
      acrophase_rad = as.numeric(phij),
      # acrophase in time units on THIS harmonic's effective period T/h
      acrophase_time = as.numeric(phij) / (2 * pi * harmonic / period),
      stringsAsFactors = FALSE),
    tests = data.frame(
      parameter = c("MESOR", "Amplitude", "Acrophase"),
      F = c(fM[1], fA[1], fP[1]), df1 = df1, df2 = df2,
      p = c(fM[2], fA[2], fP[2]), stringsAsFactors = FALSE),
    pooled = list(s2_beta = s2b, s2_gamma = s2g, cov_beta_gamma = sbg,
                  s2_mesor = s2M, s2_along = s2_along, s2_perp = s2_perp,
                  mean_direction_rad = phi_t %% (2 * pi)),
    # Bingham's own caution, as a value rather than a sentence to remember
    acrophase_differs = acro_differs,
    amplitude_interpretable = !acro_differs)
}
