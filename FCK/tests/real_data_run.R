# ==============================================================================
# tests/real_data_run.R — the audited pipeline on the real Circaflex data
#
# Deliverable 4 of the cosinor audit, now that the data and the packages are
# available. This extracts the app's OWN fitters out of server/72_harmonic.R
# (by walking the parse tree, so it is the shipped code and not a copy) and
# runs the pipeline end to end, producing the old and the new report from the
# same fits.
#
# Usage:
#   Rscript tests/real_data_run.R <path-to-xlsx> [out-prefix]
# ==============================================================================

suppressWarnings(suppressMessages({
  library(readxl); library(minpack.lm)
}))

args <- commandArgs(trailingOnly = TRUE)
xlsx <- if (length(args) >= 1) args[1] else stop("give the xlsx path")
outp <- if (length(args) >= 2) args[2] else "tests/real"

app_dir <- if (dir.exists("server")) "." else "FCK"
`%||%` <- function(a, b) if (is.null(a)) b else a
source(file.path(app_dir, "server/08_helpers_cosinor.R"))
source(file.path(app_dir, "server/09_helpers_pcanova.R"))
source(file.path(app_dir, "server/03_helpers_clock.R"))

# ---- pull the real fitters out of the shipped server file -------------------
# 72_harmonic.R is a Shiny server body: sourcing it wholesale would evaluate
# reactives. Walking the parse tree and evaluating ONLY the two function
# definitions gives us the shipped code without the session.
extract_fns <- function(file, names) {
  exprs <- parse(file, encoding = "UTF-8")
  env <- new.env(parent = globalenv())
  found <- character(0)
  walk <- function(e) {
    if (is.call(e) && length(e) >= 3 &&
        (identical(e[[1]], as.name("<-")) || identical(e[[1]], as.name("="))) &&
        is.name(e[[2]]) && as.character(e[[2]]) %in% names &&
        is.call(e[[3]]) && identical(e[[3]][[1]], as.name("function"))) {
      assign(as.character(e[[2]]), eval(e[[3]], envir = env), envir = env)
      found <<- c(found, as.character(e[[2]]))
      return(invisible(NULL))
    }
    if (is.call(e) || is.pairlist(e)) for (i in seq_along(e))
      if (!is.null(e[[i]])) try(walk(e[[i]]), silent = TRUE)
    invisible(NULL)
  }
  for (e in exprs) walk(e)
  missing <- setdiff(names, found)
  if (length(missing)) stop("could not extract: ", paste(missing, collapse = ", "))
  env
}
fenv <- extract_fns(file.path(app_dir, "server/72_harmonic.R"),
                    c("fit_cosinor", "fit_cosinor_nonlinear"))
# The fitters call the audit helpers and each other, so both must be visible
# from the environment their closures see.
for (nm in c("fck_commonality", "fck_commonality_pct", "fck_zero_amplitude_test",
             "fck_amp_se", "fck_acro_se", "fck_bingham_ci", "phi_to_hours",
             "fck_rhythm_adjusted_mean", "%||%"))
  assign(nm, get(nm, envir = globalenv()), envir = fenv)
for (nm in ls(fenv)) {
  f <- get(nm, envir = fenv)
  if (is.function(f)) { environment(f) <- fenv; assign(nm, f, envir = fenv) }
}
fit_cosinor <- get("fit_cosinor", envir = fenv)
cat("Extracted the app's own fitters from server/72_harmonic.R\n")

# ---- the data ---------------------------------------------------------------
sh_time <- "slaperigheid"; sh_cov <- "Algemene variabelen"
d <- suppressWarnings(suppressMessages(read_excel(xlsx, sheet = sh_time, .name_repair = "unique")))
a <- suppressWarnings(suppressMessages(read_excel(xlsx, sheet = sh_cov,  .name_repair = "unique")))

tc <- c("8H","9H","10H","11H","12H","14H","16H","18H","20H","21H","22H","23H","0H","2H","4H","6H")
stopifnot(all(tc %in% names(d)))
Y <- as.matrix(d[, tc]); storage.mode(Y) <- "double"
clock <- c(8,9,10,11,12,14,16,18,20,21,22,23,0,2,4,6)
tlin  <- c(8,9,10,11,12,14,16,18,20,21,22,23,24,26,28,30)
period <- 24; n_h <- 2; trend <- "exp_sat"

cat(sprintf("\nData: %d rows x %d time points. DV range [%s, %s], %.2f%% missing cells.\n",
            nrow(Y), ncol(Y), fmt1(min(Y, na.rm = TRUE)), fmt1(max(Y, na.rm = TRUE)),
            100 * mean(is.na(Y))))

# ---- THE JOIN ---------------------------------------------------------------
# The two sheets are NOT in the same order and their ID columns do not agree.
# Both joins are computed and both are reported, because choosing silently
# would decide a substantive question about the data on the analyst's behalf.
g_pos <- as.character(a$AGEcategory)                       # by row position
g_id  <- as.character(a$AGEcategory)[match(d$ID, a$ID)]    # by ID
cat("\n--- THE GROUP JOIN ---\n")
cat("Row-position join: "); print(table(g_pos, useNA = "ifany"))
cat("ID join:           "); print(table(g_id,  useNA = "ifany"))
cat(sprintf("The two agree for %d of %d rows (%.1f%%).\n",
            sum(g_pos == g_id, na.rm = TRUE), nrow(d),
            100 * mean(g_pos == g_id, na.rm = TRUE)))
cat(sprintf("slaperigheid has %d rows but only %d distinct IDs (%d repeats).\n",
            nrow(d), length(unique(d$ID)), sum(duplicated(d$ID))))
cat("The app joins by ROW POSITION, which is what reproduces the reported\n")
cat("654/410/181/59. Whether that is CORRECT is a question about the file.\n")
grp <- g_pos

# ---- fit every subject with the app's own fitter -----------------------------
cat("\nFitting", nrow(Y), "subjects (exp_sat, 2 harmonics, raw data)...\n")
fits <- vector("list", nrow(Y))
t0 <- Sys.time()
for (i in seq_len(nrow(Y))) {
  fits[i] <- list(tryCatch(fit_cosinor(tlin, as.numeric(Y[i, ]), period = period,
                                       n_harmonics = n_h, trend_type = trend),
                           error = function(e) list(success = FALSE, message = conditionMessage(e))))
  if (i %% 200 == 0) cat("  ", i, " (", round(as.numeric(Sys.time() - t0, units = "secs")), "s)\n", sep = "")
}
okv  <- vapply(fits, function(f) isTRUE(f$success), logical(1))
conv <- vapply(fits, function(f) isTRUE(f$success) && isTRUE(f$converged), logical(1))
bnd  <- vapply(fits, function(f) isTRUE(f$success) && isTRUE(f$boundary_hit), logical(1))
cat(sprintf("\nreturned a fit: %d   converged: %d   on a bound: %d   failed: %d\n",
            sum(okv), sum(conv & !bnd), sum(bnd), sum(!okv)))

keep_old <- which(okv)              # what the old code summarised
keep_new <- which(conv & !bnd)      # what the audited code summarises

P <- function(idx) {
  f <- fits[idx]
  data.frame(
    row = idx,
    intercept = vapply(f, function(z) as.numeric(z$mesor), 1),
    A_sat = vapply(f, function(z) as.numeric(z$trend_params$A_sat$coef %||% NA), 1),
    tau   = vapply(f, function(z) as.numeric(z$trend_params$tau$coef %||% NA), 1),
    amp1  = vapply(f, function(z) z$amplitudes[1], 1),
    phi1  = vapply(f, function(z) z$acrophases[1], 1),
    amp2  = vapply(f, function(z) z$amplitudes[2], 1),
    phi2  = vapply(f, function(z) z$acrophases[2], 1),
    acro2_time = vapply(f, function(z) z$acrophases_time[2], 1),
    r2    = vapply(f, function(z) z$r_squared, 1),
    r2_S  = vapply(f, function(z) z$r_squared_S, 1),
    r2_C  = vapply(f, function(z) z$r_squared_C, 1),
    uS    = vapply(f, function(z) z$unique_S %||% NA_real_, 1),
    uC    = vapply(f, function(z) z$unique_C %||% NA_real_, 1),
    sh    = vapply(f, function(z) z$shared_SC %||% NA_real_, 1),
    p     = vapply(f, function(z) z$p_value, 1),
    grp   = grp[idx], stringsAsFactors = FALSE)
}
Po <- P(keep_old); Pn <- P(keep_new)

summ <- function(D) {
  r1 <- fck_resultants(D$phi1, D$amp1); r2 <- fck_resultants(D$phi2, D$amp2)
  list(D = D, r1 = r1, r2 = r2,
       int = mean(D$intercept, na.rm = TRUE),
       A = mean(D$A_sat, na.rm = TRUE), tau = mean(D$tau, na.rm = TRUE),
       ram = fck_rhythm_adjusted_mean(mean(D$intercept, na.rm = TRUE), trend,
                                      c(mean(D$A_sat, na.rm = TRUE), mean(D$tau, na.rm = TRUE)),
                                      min(tlin), max(tlin), min(tlin)))
}
So <- summ(Po); Sn <- summ(Pn)

line <- function() cat(strrep("=", 78), "\n")

# =============================== OLD =========================================
sink(paste0(outp, "_report_OLD.txt"))
line(); cat("CIRCAFLEX 2019 — REPORT AS THE CODE PRODUCED IT BEFORE THE AUDIT\n"); line()
cat("\nPeriod: 24   Harmonics: 2   Trend: SATURATING EXPONENTIAL\n")
cat("Number of subjects:", nrow(Y), "\n")
cat("  - Successfully fitted:", length(keep_old), "\n")
cat("\n--- Population Mean Parameters (Vector-Averaged) ---\n")
cat(sprintf("MESOR:     %.3f\n", So$int))
cat(sprintf("Amplitude (H1): %.3f\n", So$r1$mean_amplitude))
cat(sprintf("Acrophase (H1): %.2f (%.2f hours)\n",
            So$r1$mean_dir_weighted * 180 / pi, So$r1$mean_dir_weighted * period / (2 * pi)))
cat("\n  All Harmonics (vector-averaged):\n")
cat(sprintf("    H1: Amplitude=%.3f, Acrophase=%.2f hours\n",
            So$r1$mean_amplitude, So$r1$mean_dir_weighted * period / (2 * pi)))
cat(sprintf("    H2: Amplitude=%.3f, Acrophase=%.2f hours\n",
            So$r2$mean_amplitude, So$r2$mean_dir_weighted * period / (2 * pi) / 2))
z_old <- nrow(Po) * So$r1$r_weighted^2
cat("\nRayleigh test for uniformity (H1):\n")
cat(sprintf("  Z = %.3f, p = %.4f\n", z_old, exp(-z_old)))
cat(sprintf("  Mean resultant length (r-bar) = %.3f\n", So$r1$r_weighted))
cat("\n--- Fitted Model Equation ---\n")
cat(sprintf("Y(t) = %.2f + %.2f*cos(2pi*1*t/24 - %.2f) + %.2f*cos(2pi*2*t/24 - %.2f)\n",
            So$int, So$r1$mean_amplitude, So$r1$mean_dir_weighted,
            So$r2$mean_amplitude, So$r2$mean_dir_weighted))
cat("\n--- Individual Parameter Summary ---\n")
cat(sprintf("MESOR:     Mean=%.3f, SD=%.3f\n", mean(Po$intercept), sd(Po$intercept)))
cat(sprintf("A_sat (asymptote): Mean=%.3f, SD=%.3f (units)\n",
            mean(Po$A_sat, na.rm = TRUE), sd(Po$A_sat, na.rm = TRUE)))
cat(sprintf("tau (time constant): Mean=%.2f, SD=%.2f (hours)\n",
            mean(Po$tau, na.rm = TRUE), sd(Po$tau, na.rm = TRUE)))
cat(sprintf("H1 Amplitude: Mean=%.3f, SD=%.3f\n", mean(Po$amp1), sd(Po$amp1)))
cat(sprintf("H1 Acrophase: Circular mean=%.2f h, r-bar=%.3f\n",
            So$r1$mean_dir_unweighted * period / (2 * pi), So$r1$r_unweighted))
cat(sprintf("H2 Acrophase: Arithmetic mean=%.2f h   [on a 0-24 scale: the /h divisor was missing]\n",
            mean(Po$acro2_time * 2, na.rm = TRUE)))
cat(sprintf("R-squared: Mean=%.3f, Range=[%.3f, %.3f]\n",
            mean(Po$r2), min(Po$r2), max(Po$r2)))
cat(sprintf("Significant rhythms (p<0.05): %d / %d (%.1f%%)\n",
            sum(Po$p < 0.05, na.rm = TRUE), nrow(Po),
            100 * mean(Po$p < 0.05, na.rm = TRUE)))
cat("\n--- Variance Decomposition: Relative Importance ---\n")
ps <- mean(100 * Po$r2_S / Po$r2, na.rm = TRUE); pc <- mean(100 * Po$r2_C / Po$r2, na.rm = TRUE)
cat(sprintf("R2 from Process S (homeostatic): Mean=%.3f\n", mean(Po$r2_S, na.rm = TRUE)))
cat(sprintf("R2 from Process C (circadian):   Mean=%.3f\n", mean(Po$r2_C, na.rm = TRUE)))
cat(sprintf("\nProportion of total R2 explained by:\n  Process S: %.1f%%\n  Process C: %.1f%%\n", ps, pc))
cat(sprintf("  [these sum to %.1f%%]\n", ps + pc))
cat(sprintf("\nInterpretation: %s is the dominant component\n",
            if (pc > ps) "Circadian rhythm (C)" else "Homeostatic process (S)"))
cat("\n--- Group-Specific Parameters ---\n")
gs <- 0L
for (g in c("YOUTH","ADULT","MIDDLE_AGE","ELDERLY")) {
  idx <- which(!is.na(Po$grp) & Po$grp == g); if (length(idx) < 3) next
  gs <- gs + length(idx)
  rr <- fck_resultants(Po$phi1[idx], Po$amp1[idx])
  cat(sprintf("\nGroup '%s' (n=%d):\n", g, length(idx)))
  cat(sprintf("  MESOR:     %.3f (SD=%.3f)\n", mean(Po$intercept[idx]), sd(Po$intercept[idx])))
  cat(sprintf("  H1 Amplitude: %.3f (SD=%.3f)\n", rr$mean_amplitude, sd(Po$amp1[idx])))
  cat(sprintf("  H1 Acrophase: %.2f hours\n", rr$mean_dir_weighted * period / (2 * pi)))
}
cat(sprintf("\n[group n's sum to %d; %d subjects were fitted]\n", gs, nrow(Po)))
sink()

# =============================== NEW =========================================
sink(paste0(outp, "_report_NEW.txt"))
line(); cat("CIRCAFLEX 2019 — REPORT AFTER THE AUDIT\n"); line()
cat("\nDependent variable: sleepiness VAS\n")
cat(sprintf("  Observed range in the file: [%s, %s]\n", fmt1(min(Y, na.rm=TRUE)), fmt1(max(Y, na.rm=TRUE))))
cat("  Admissible range: NOT SPECIFIED. Set it in the DV panel to have fitted\n")
cat("  values checked; on a 0-500 VAS a fit predicting -20 is misspecified, and\n")
cat("  nothing in the report can tell you that until the bounds are declared.\n")
cat(sprintf("Period: 24 h   Harmonics: 2   Trend: SATURATING EXPONENTIAL\n"))
cat("Data: RAW (no smoothing applied)\n")
cat(sprintf("  %.2f%% of cells are missing; the cosinor uses the observations present.\n",
            100 * mean(is.na(Y))))
cat("Time origin: MIDNIGHT. The trend is anchored at the first observation (t - 8)\n")
cat("  while the harmonics are anchored at midnight, so the intercept is the\n")
cat("  constant of a model with two origins and is not the value at either.\n")

cat("\n--- Fit outcomes ---\n")
cat(sprintf("Subjects attempted:         %d\n", nrow(Y)))
cat(sprintf("  Converged:                %d\n", sum(conv & !bnd)))
cat(sprintf("  Converged but on a bound: %d  (excluded)\n", sum(bnd)))
cat(sprintf("  Failed outright:          %d\n", sum(!okv)))
cat(sprintf("Population summaries below use %d subject(s).\n", nrow(Pn)))
cat(sprintf("The old report counted all %d as successes.\n", sum(okv)))

cat("\n--- Central value ---\n")
cat(sprintf("Intercept (beta_0, at t = 0):                           %s\n", fmt3(Sn$int)))
cat(sprintf("MESOR (rhythm-adjusted mean over the observed window):  %s\n", fmt3(Sn$ram)))
cat(sprintf("  Difference: %s units. The intercept is NOT the MESOR here.\n", fmt2(Sn$ram - Sn$int)))

cat("\n--- Population rhythm parameters (VECTOR-averaged) ---\n")
cat(sprintf("  H1: amplitude = %s, acrophase = %s h\n",
            fmt3(Sn$r1$mean_amplitude), fmt2(phi_to_hours(Sn$r1$mean_dir_weighted, period, 1))))
cat(sprintf("  H2: amplitude = %s, acrophase = %s h  (modulo 12 h)\n",
            fmt3(Sn$r2$mean_amplitude), fmt2(phi_to_hours(Sn$r2$mean_dir_weighted, period, 2))))

cat("\n--- Circular concentration and the Rayleigh test ---\n")
for (h in 1:2) {
  rr <- if (h == 1) Sn$r1 else Sn$r2
  ry <- fck_rayleigh(rr$r_unweighted, rr$n)
  cat(sprintf("H%d  n = %d\n", h, rr$n))
  cat(sprintf("  r-bar (UNWEIGHTED, for Rayleigh):                %s\n", fmt3(rr$r_unweighted)))
  cat(sprintf("  r-bar (AMPLITUDE-WEIGHTED, for the vector mean): %s\n", fmt3(rr$r_weighted)))
  cat(sprintf("  Rayleigh Z = %s, p = %s\n", fmt1e(ry$Z),
              format.pval(ry$p, digits = 3, eps = 1e-16)))
  cat(sprintf("  Circular SD = %s h\n", fmt3(phi_to_hours(rr$circ_sd_rad, period, h))))
}

cat("\n--- Fitted model equation (pooled) ---\n")
cat(fck_format_equation(Sn$int, trend, c(Sn$A, Sn$tau),
                        c(Sn$r1$mean_amplitude, Sn$r2$mean_amplitude),
                        c(Sn$r1$mean_dir_weighted, Sn$r2$mean_dir_weighted),
                        period, min(tlin)), "\n")

cat("\n--- Individual parameters (ARITHMETIC means, +/- linear SD) ---\n")
cat(sprintf("  Intercept:           %s (SD %s)   [arithmetic]\n",
            fmt3(mean(Pn$intercept)), fmt3(sd(Pn$intercept))))
cat(sprintf("  A_sat (asymptote):   %s (SD %s)   [arithmetic]\n",
            fmt3(mean(Pn$A_sat, na.rm=TRUE)), fmt3(sd(Pn$A_sat, na.rm=TRUE))))
cat(sprintf("  tau (time constant): %s (SD %s) h   [arithmetic]\n",
            fmt2(mean(Pn$tau, na.rm=TRUE)), fmt2(sd(Pn$tau, na.rm=TRUE))))
if (sd(Pn$tau, na.rm=TRUE) / mean(Pn$tau, na.rm=TRUE) > 0.5)
  cat("    ! SD/mean above 0.5: a likelihood ridge, not population heterogeneity.\n")
for (h in 1:2) {
  rr <- if (h == 1) Sn$r1 else Sn$r2
  ac <- if (h == 1) Pn$amp1 else Pn$amp2
  cat(sprintf("  H%d amplitude: %s (SD %s)   [arithmetic]\n", h, fmt3(mean(ac)), fmt3(sd(ac))))
  cat(sprintf("  H%d acrophase: circular mean %s h, circular SD %s h   [circular]\n",
              h, fmt2(phi_to_hours(rr$mean_dir_unweighted, period, h)),
              fmt2(phi_to_hours(rr$circ_sd_rad, period, h))))
}
cat(sprintf("\n  R-squared: mean %s, range [%s, %s]\n",
            fmt3(mean(Pn$r2)), fmt3(min(Pn$r2)), fmt3(max(Pn$r2))))
cat(sprintf("  Significant rhythms (p<0.05): %d / %d (%s%%)\n",
            sum(Pn$p < 0.05, na.rm = TRUE), nrow(Pn), fmt1(100 * mean(Pn$p < 0.05, na.rm = TRUE))))
cat("    The zero-amplitude test is now the harmonics GIVEN the trend.\n")

cat("\n--- Variance decomposition (commonality analysis) ---\n")
uS <- mean(Pn$uS, na.rm=TRUE); uC <- mean(Pn$uC, na.rm=TRUE); shd <- mean(Pn$sh, na.rm=TRUE)
tot <- mean(Pn$r2)
cat(sprintf("  Unique to Process S (homeostatic): %s  (%s%%)\n", fmt3(uS), fmt1(100*uS/tot)))
cat(sprintf("  Unique to Process C (circadian):   %s  (%s%%)\n", fmt3(uC), fmt1(100*uC/tot)))
cat(sprintf("  Shared between S and C:            %s  (%s%%)\n", fmt3(shd), fmt1(100*shd/tot)))
cat(          "  ------------------------------------------------\n")
cat(sprintf("  Total R-squared:                   %s  (%s%%)\n",
            fmt3(uS+uC+shd), fmt1(100*(uS+uC+shd)/tot)))
cat("  No dominance verdict is emitted.\n")

cat("\n--- Group-specific parameters ---\n")
ga <- fck_group_audit(Pn$grp, seq_len(nrow(Pn)))
cat(sprintf("Group sizes sum to %d; %d subjects were fitted.\n", sum(ga$counts), nrow(Pn)))
if (ga$n_unassigned > 0)
  cat(sprintf("  %d subject(s) have no usable group label and appear as UNASSIGNED below.\n",
              ga$n_unassigned))
cat("  These reconcile.\n\n")
for (g in c("YOUTH","ADULT","MIDDLE_AGE","ELDERLY","__UNASSIGNED__")) {
  idx <- if (g == "__UNASSIGNED__") which(is.na(Pn$grp)) else which(!is.na(Pn$grp) & Pn$grp == g)
  if (!length(idx)) next
  gI <- mean(Pn$intercept[idx]); gA <- mean(Pn$A_sat[idx], na.rm=TRUE); gT <- mean(Pn$tau[idx], na.rm=TRUE)
  gram <- fck_rhythm_adjusted_mean(gI, trend, c(gA, gT), min(tlin), max(tlin), min(tlin))
  rr <- fck_resultants(Pn$phi1[idx], Pn$amp1[idx]); r2g <- fck_resultants(Pn$phi2[idx], Pn$amp2[idx])
  cat(sprintf("Group '%s' (n = %d):\n",
              if (g == "__UNASSIGNED__") "UNASSIGNED (no usable group label)" else g, length(idx)))
  if (g == "__UNASSIGNED__")
    cat("  Shown because these subjects exist. They are not a group; do not compare them.\n")
  cat(sprintf("  Intercept (beta_0, at t = 0):  %s (SD %s)   [arithmetic]\n",
              fmt3(gI), fmt3(sd(Pn$intercept[idx]))))
  cat(sprintf("  MESOR (rhythm-adjusted mean):  %s   [integrated over the window]\n", fmt3(gram)))
  cat(sprintf("  A_sat: %s   tau: %s h   [arithmetic]\n", fmt3(gA), fmt2(gT)))
  cat(sprintf("  H1 amplitude: %s   [vector]   /  %s (SD %s)   [arithmetic]\n",
              fmt3(rr$mean_amplitude), fmt3(mean(Pn$amp1[idx])), fmt3(sd(Pn$amp1[idx]))))
  cat(sprintf("  H1 acrophase: %s h   [vector]   r-bar unweighted %s, weighted %s\n",
              fmt2(phi_to_hours(rr$mean_dir_weighted, period, 1)),
              fmt3(rr$r_unweighted), fmt3(rr$r_weighted)))
  cat("\n  Fitted equation:\n  ")
  cat(fck_format_equation(gI, trend, c(gA, gT),
                          c(rr$mean_amplitude, r2g$mean_amplitude),
                          c(rr$mean_dir_weighted, r2g$mean_dir_weighted),
                          period, min(tlin)), "\n\n")
}

cat("--- Group comparisons ---\n")
lf <- factor(Pn$grp, levels = c("YOUTH","ADULT","MIDDLE_AGE","ELDERLY"))
for (nm in c("intercept","amp1","amp2","A_sat","tau")) {
  lt <- fck_group_linear_test(Pn[[nm]], lf); tr <- fck_group_trend_test(Pn[[nm]], lf)
  if (is.null(lt)) next
  cat(sprintf("\n%s:\n", nm))
  cat(sprintf("  F(%d, %d) = %s, p = %s;  eta^2 = %s, omega^2 = %s\n",
              lt$df1, lt$df2, fmt3(lt$F), format.pval(lt$p, digits=3, eps=1e-16),
              fmt3(lt$eta2), fmt3(lt$omega2)))
  for (gg in lt$levels)
    cat(sprintf("    %-12s %s (SD %s, n = %d)\n", gg, fmt3(lt$means[gg]), fmt3(lt$sds[gg]), lt$ns[gg]))
  if (!is.null(lt$largest))
    cat(sprintf("  Largest contrast %s - %s: %s, 95%% CI [%s, %s], Hedges' g = %s\n",
                lt$largest$a, lt$largest$b, fmt3(lt$largest$diff),
                fmt3(lt$largest$ci[1]), fmt3(lt$largest$ci[2]), fmt3(lt$largest$hedges_g)))
  if (!is.null(tr))
    cat(sprintf("  Monotone trend (YOUTH<ADULT<MIDDLE_AGE<ELDERLY): t(%d) = %s, p = %s, 95%% CI [%s, %s]\n",
                tr$df, fmt3(tr$t), format.pval(tr$p, digits=3, eps=1e-16),
                fmt3(tr$ci[1]), fmt3(tr$ci[2])))
}
sink()

# =============================== the diff ====================================
line(); cat("OLD vs NEW, on the real data\n"); line()
row <- function(what, old, new) cat(sprintf("%-46s %14s %14s\n", what, old, new))
cat(sprintf("%-46s %14s %14s\n", "", "OLD", "NEW"))
row("subjects in the population summary", length(keep_old), nrow(Pn))
row("central value reported as the level", fmt2(So$int), fmt2(Sn$ram))
row("Rayleigh Z (H1)", fmt1e(nrow(Po) * So$r1$r_weighted^2),
    fmt1e(fck_rayleigh(Sn$r1$r_unweighted, Sn$r1$n)$Z))
row("r-bar used for Rayleigh", fmt3(So$r1$r_weighted), fmt3(Sn$r1$r_unweighted))
row("H2 acrophase (h)", fmt2(mean(Po$acro2_time * 2, na.rm = TRUE)),
    fmt2(phi_to_hours(Sn$r2$mean_dir_weighted, period, 2)))
row("variance parts sum to (% of total R2)",
    fmt1(mean(100*Po$r2_S/Po$r2, na.rm=TRUE) + mean(100*Po$r2_C/Po$r2, na.rm=TRUE)),
    fmt1(100*(uS+uC+shd)/tot))
row("group n's sum to", sum(table(Po$grp[!is.na(Po$grp)])), sum(ga$counts) + ga$n_unassigned)
row("significant rhythms (%)",
    fmt1(100*mean(Po$p < 0.05, na.rm=TRUE)), fmt1(100*mean(Pn$p < 0.05, na.rm=TRUE)))

y27 <- function(with_trend) {
  v <- Sn$int + Sn$r1$mean_amplitude * cos(2*pi*27/period - Sn$r1$mean_dir_weighted) +
    Sn$r2$mean_amplitude * cos(2*pi*2*27/period - Sn$r2$mean_dir_weighted)
  if (with_trend) v <- v + Sn$A * (1 - exp(-(27 - min(tlin)) / Sn$tau))
  v
}
obs27 <- mean(Y[, which(tlin == 26)], na.rm = TRUE)
cat(sprintf("\nPooled equation at t = 27 h (03:00 next day):\n"))
cat(sprintf("  OLD printed equation:  %s\n", fmt2(y27(FALSE))))
cat(sprintf("  NEW printed equation:  %s\n", fmt2(y27(TRUE))))
cat(sprintf("  Observed mean at 02:00 (nearest measured point): %s\n", fmt2(obs27)))
cat(sprintf("  The old equation under-predicted by %s units.\n", fmt2(y27(TRUE) - y27(FALSE))))

saveRDS(list(fits_summary = Pn, grp_pos = g_pos, grp_id = g_id, Y = Y, tlin = tlin),
        paste0(outp, "_fits.rds"))
line()
cat("Wrote ", outp, "_report_OLD.txt and ", outp, "_report_NEW.txt\n", sep = "")
