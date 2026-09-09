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

# ==============================================================================
# THE JOINT TEST ON THE RHYTHMIC VECTOR (P20/R4)
# ==============================================================================
# One-way MANOVA on the bivariate (cosine, sine) coefficient pair. This is the
# test the marginal amplitude and acrophase F ratios are marginals OF, and it is
# the only one of the three that can see an antipodal difference (see the note
# on the blind spot below). It lives here, as a pure function, because
# server/72_harmonic.R needs exactly the same test for Bingham's parameter test
# and a second implementation of a standard procedure is how the two drift
# apart.
#
# For p = 2 response variables Rao's transformation of Wilks' Lambda is EXACT,
# and reduces to the two-sample Hotelling T-squared when there are two groups.
# stats::manova is R's own implementation of it, so the arithmetic is not
# re-derived here.
dance_bivariate_manova <- function(x, y, group) {
  ok <- is.finite(x) & is.finite(y) & !is.na(group)
  x <- x[ok]; y <- y[ok]; group <- droplevels(factor(group[ok]))
  gl <- levels(group); k <- length(gl); N <- length(x)
  if (k < 2) return(list(ok = FALSE, message = "Need at least 2 groups."))
  ns <- as.integer(table(group))
  if (any(ns < 3)) return(list(ok = FALSE,
    message = "Each group needs at least 3 participants for the joint vector test."))
  Y <- cbind(x, y)
  E <- Reduce("+", lapply(gl, function(g) {
    m <- Y[group == g, , drop = FALSE]
    cm <- colMeans(m)
    t(sweep(m, 2, cm)) %*% sweep(m, 2, cm)
  }))
  if (!is.finite(det(E)) || det(E) < .Machine$double.eps)
    return(list(ok = FALSE, message = "Singular within-group covariance."))
  fit <- tryCatch(summary(stats::manova(Y ~ group), test = "Wilks"),
                  error = function(e) NULL)
  if (is.null(fit) || is.null(fit$stats))
    return(list(ok = FALSE, message = "The joint vector test could not be fitted."))
  st <- fit$stats["group", ]
  list(ok = TRUE, F = unname(st[["approx F"]]), df1 = unname(st[["num Df"]]),
       df2 = unname(st[["den Df"]]), p = unname(st[["Pr(>F)"]]),
       lambda = unname(st[["Wilks"]]), n_groups = k, n_total = N)
}

# Largest pairwise circular separation among a set of angles, in [0, pi].
# Used to decide whether the LINEARISED acrophase test below has any power.
dance_max_angular_sep <- function(phi) {
  phi <- phi[is.finite(phi)]
  if (length(phi) < 2) return(0)
  d <- outer(phi, phi, "-")
  max(abs(atan2(sin(d), cos(d))))
}

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

  # ---- the blind spot of the linearised acrophase test (P20/R4) -------------
  # sin^2(phi_j - phi_tilde) is the squared PERPENDICULAR displacement of group
  # j from the pooled direction, and it is pi-periodic: two groups half a cycle
  # apart lie on the same line through the origin, so every displacement is
  # zero. With group acrophases at 0 h and 12 h and n = 20 each, the test
  # returned F = 0.011, p = .918 -- and, because the amplitude test compares the
  # UNSIGNED lengths |A_j|, which are also equal, that marginal was null too.
  # The largest acrophase difference the design can express was invisible to
  # both marginals, and `amplitude_interpretable` then said TRUE, which is the
  # exact opposite of Bingham's caution.
  #
  # WHAT IS AND IS NOT SALVAGED. The failure is one of POWER, not of level: a
  # large angular separation SHRINKS the numerator, so the test stays valid --
  # a significant acrophase result is a real rejection at any separation. What
  # is not supportable is reading a NULL result as "the acrophases agree". The
  # statistic is monotone in the separation only up to a quarter cycle; past
  # pi/2 it turns around and heads back to zero. So pi/2 is where the claim is
  # cut off, and the flag below says so rather than leaving it implicit.
  #
  # The gap is covered rather than merely declared: the joint MANOVA on the
  # (beta, gamma) pair -- of which amplitude and acrophase are marginals -- has
  # no such blind spot, and is returned alongside. It is the test to read when
  # the linearisation is unsupported.
  sep <- dance_max_angular_sep(as.numeric(phij))
  acro_supported <- is.finite(sep) && sep <= pi / 2
  joint <- dance_bivariate_manova(beta, gamma, group)

  acro_differs <- is.finite(fP[2]) && fP[2] < 0.05
  # A joint difference that the marginals cannot resolve still means the groups'
  # rhythmic vectors differ, and an amplitude comparison across vectors that
  # point different ways is what Bingham warns against. So the amplitude is
  # declared interpretable only when the acrophase claim is actually supported.
  joint_differs <- isTRUE(joint$ok) && is.finite(joint$p) && joint$p < 0.05
  amp_ok <- (!acro_differs) && acro_supported && !(joint_differs && !acro_differs && sep > pi / 4)

  acro_note <- if (acro_supported) NULL else sprintf(paste0(
    "The group acrophases are up to %.1f h apart (%.2f rad) on this harmonic's ",
    "%.1f h effective period. Bingham's acrophase F ratio measures displacement ",
    "PERPENDICULAR to the pooled direction, which is half-cycle periodic: past a ",
    "quarter cycle it decreases again and at half a cycle it is exactly zero. A ",
    "significant result would still be a real difference, but a NULL result here ",
    "carries no information. Read the joint test on the (cosine, sine) vector ",
    "instead, and compare the per-group acrophases descriptively."),
    sep / (2 * pi * harmonic / period), sep, period / harmonic)

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
    # the joint test on the vector: no blind spot, and the one to read when the
    # linearised marginals are unsupported
    joint = joint,
    max_angular_sep_rad = sep,
    max_angular_sep_time = sep / (2 * pi * harmonic / period),
    acrophase_test_supported = acro_supported,
    acrophase_test_note = acro_note,
    # Bingham's own caution, as a value rather than a sentence to remember
    acrophase_differs = acro_differs,
    amplitude_interpretable = amp_ok)
}
