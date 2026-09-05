# ==============================================================================
# server/01b_kernel_contracts.R — the statistical contracts, in one place
# ==============================================================================
# AUDIT (P2.3). Both reviews made the same point: the app has unusually good
# audit comments, but the statistical CONTRACTS live inside 5,800-line server
# files, so the only way to learn what a kernel assumes is to read its body.
#
# Extracting the kernels into a package is the right long-term move and is not
# this change: moving eight estimators out of the reactive files at once is how
# a working app stops working. What is recorded here is the contract for each
# one — orientation, units, estimand, assumptions, missing-data policy, what
# the intervals mean, and the failure modes — as documentation that sits beside
# the code and is checked by tests/testthat/test-p2-corrections.R.
#
# The contracts are deliberately specific about the things that were WRONG at
# some point, because those are the ones a future reader will get wrong again.
#
# Format: FCK_CONTRACTS[[name]] is a list. `orientation` is stated for every
# kernel because the app mixes both conventions and that is the single easiest
# way to introduce a silent transpose.
# ==============================================================================

FCK_CONTRACTS <- list(

  fit_cosinor = list(
    file        = "server/72_harmonic.R",
    model       = "y(t) = MESOR + trend(t) + sum_h [b_cos_h cos(w_h t) + b_sin_h sin(w_h t)] + e",
    orientation = "one subject at a time: y and time are vectors of equal length",
    time_units  = "MODEL time (hours since the chosen origin), never clock time. Clock conversion happens only for display, through fck_model_to_clock().",
    estimand    = "per-subject MESOR, trend parameters, and amplitude/acrophase per harmonic",
    estimator   = "OLS on the linear-in-parameters design",
    assumptions = "Gaussian errors, independent within subject, trend correctly specified",
    missing     = "rows with NA in y or time are dropped before fitting",
    bounds      = "none; bounds apply to the nonlinear fitter only",
    intervals   = "parametric, from the OLS covariance; amplitude and acrophase SEs by the delta method",
    zero_amp    = "F test of the harmonics GIVEN the trend: full model against a refit trend-only model, df1 = 2*n_harmonics",
    degenerate  = "constant y gives SST = 0, so R^2 is NA (fck_r_squared), not NaN; a perfect fit floors sigma^2 (fck_gaussian_loglik) so AIC stays finite",
    returns     = "list(success, mesor, amplitudes, acrophases, acrophases_time, r_squared, p_value, aic, ...)"
  ),

  fit_cosinor_nonlinear = list(
    file        = "server/72_harmonic.R",
    model       = "as fit_cosinor, with trend = A_sat * (1 - exp(-(t - t0)/tau)) for exp_sat",
    orientation = "one subject at a time",
    time_units  = "model time; t_offset is the first observation, not midnight",
    estimand    = "as fit_cosinor, plus A_sat and tau",
    estimator   = "minpack.lm::nlsLM, falling back to nls(algorithm='port'); convergence is READ from the fit object, never assumed from the absence of an error",
    assumptions = "Gaussian errors; the exp_sat rise is monotone and saturating by design (extended wakefulness: there is no sleep episode in the window, so a piecewise rise/decay is NOT the right specification here)",
    bounds      = "amplitude_max constrains sqrt(b_cos^2 + b_sin^2), enforced by refitting inside the inscribed box only when a converged fit violates it; amplitude_min is REPORTED, never imposed, because forcing a floor would invent a rhythm",
    zero_amp    = "the reduced model REFITS tau (fck_reduced_exp_sat_sse). Freezing tau at its full-model value inflated Type-I error from 5% to 14.6%",
    boundary    = "tau frequently hits its bound on real data; asymptotic F and CIs are questionable there and the bounds table reports which fits are affected",
    degenerate  = "as fit_cosinor",
    returns     = "as fit_cosinor, plus trend_params, amp_bound_action, amplitude_below_min, convergence"
  ),

  fck_auto_lambda = list(
    file        = "server/04_helpers_fd.R",
    estimand    = "one smoothing parameter for the whole sample",
    orientation = "data_mat is subjects x time",
    estimator   = "coarse log-grid then optimize() on mean GCV, using fda::smooth.basis -- the SAME estimator that then does the smoothing",
    rationale   = "one lambda for the sample, not per subject: a per-subject lambda makes curves incomparable, which is the premise of everything downstream in the fPCA family",
    assumptions = "GCV is a reasonable criterion for these data; the grid brackets the optimum",
    missing     = "subjects with fewer than min_points_needed observed values are excluded from the objective",
    failure     = "returns NULL when no subject can be scored; the caller must say so rather than silently using lambda = 0",
    warning     = "lambda = 0 is NOT automatic selection. fda::smooth.basis does not search: it smooths with the lambda given and reports that lambda's GCV. With nbasis = n_time, lambda = 0 interpolates."
  ),

  perform_functional_anova = list(
    file        = "server/50_fanova.R",
    orientation = "curves is TIME x CURVES (n_time rows). The RM kernel builds subjects x visits matrices per time point -- the two conventions differ and a transpose here is silent.",
    estimand    = "equality of group mean functions",
    estimator   = "pointwise one-way F at each evaluation point",
    inference   = "permutation of group labels; p = (1 + #{F* >= F}) / (1 + B) so a Monte Carlo p is never 0; FDR-adjusted across time points",
    exchange    = "group labels are exchangeable under the null: between-subject design, whole curves permuted",
    global      = "L2 norm of SSB/n, integrated by trapezoid over the evaluation grid (fck_l2_norm) -- NOT sum of squares over grid points, which scales with grid density",
    intervals   = "POINTWISE bootstrap percentile intervals, 2000 replicates. Not simultaneous bands: a region where the interval excludes zero is not a family-wise claim.",
    missing     = "complete curves only",
    returns     = "list(F_stat, p_values_pointwise, p_values_adjusted, L2_stat, p_value_L2, sig_regions, ...)"
  ),

  perform_rm_fanova = list(
    file        = "server/50_fanova.R",
    orientation = "internally builds one subjects x visits matrix PER TIME POINT; the complete-case set can differ between time points, so Y_rows carries the subject indices behind each",
    estimand    = "equality of condition mean functions within subject",
    estimator   = "pointwise one-factor repeated-measures F. SS_error removes BOTH margins and adds the grand mean back; df = (n-1)(k-1). Removing only the subject margin leaves the condition effect in the residual and returns 22.8% of the correct F.",
    inference   = "within-subject permutation of the condition labels, ONE relabelling per subject per replicate applied across the whole trajectory",
    exchange    = "condition labels are exchangeable within a subject under the null; the exchangeable unit is the whole curve, not the value at a time point",
    global      = "none. p_value_L2 is NA: no valid global statistic is computed, and the pointwise p-values are what the app reports.",
    package     = "this does NOT call rmfanova. The app previously guessed four signatures for it, none of which is that package's API (it takes a list of condition-specific n x p matrices and requires a complete balanced design).",
    missing     = "complete cases per time point",
    returns     = "list(F_stat, p_values_pointwise, p_values_adjusted, df_between, df_within, ...)"
  ),

  linear_shift_alignment = list(
    file        = "server/40_fpca.R",
    orientation = "curves is TIME x CURVES",
    estimand    = "one time shift per curve, against the sample mean (or median, or first curve)",
    estimator   = "cross-correlation lag, used AS MEASURED. Two attenuation constants (0.1 then 0.5) previously shrank it to 5% of itself.",
    warp        = "h(t) = t - s, monotone by construction, deterministic. No random term: a runif() perturbation was previously added to the estimate under a comment about visualisation.",
    limits      = "|s| is capped at a quarter of the domain; beyond that a shift is not identified",
    determinism = "same input, same output, bit for bit. There is no set.seed in this module because there is no longer anything to seed.",
    validation  = "tests/warping_test.R recovers a known phase shift at r = 0.998, slope 0.901; identity in gives identity out exactly",
    caveat      = "this is a SHIFT registration. It cannot represent differential stretching. For that, use a monotone-spline or SRVF/Fisher-Rao method; the parametric family here is restricted to alpha ranges where the map is a bijection."
  ),

  fit_fosr = list(
    file        = "server/70_fosr.R",
    orientation = "Y is SUBJECTS x TIME; coefficient curves are terms x time",
    estimand    = "beta(t), the effect of each scalar predictor on the response curve",
    estimator   = "pointwise OLS via QR with an explicit rank check; the normal equations are not inverted",
    inference   = "studentised statistic, FDR-adjusted across time points",
    bootstrap   = "residual bootstrap, used for the SE and the percentile interval ONLY. It resamples around the fitted ALTERNATIVE, so it is not a null distribution: 2*mean(boot <= 0) is CI inversion, not a hypothesis test, and is no longer reported as a p-value.",
    factors     = "in the GAM branch, a factor emits one contrast curve per non-reference level; they were previously left as rows of zeros",
    missing     = "complete cases",
    failure     = "a rank-deficient design is an error naming the collinear columns, not a silent pseudo-inverse"
  ),

  fit_sofr = list(
    file        = "server/71_sofr.R",
    estimand    = "a scalar outcome regressed on a predictor CURVE, via refund::pfr",
    orientation = "X_func is SUBJECTS x TIME",
    argvals     = "the real normalised time positions are passed to lf(), so beta(t) is a function of time and not of column index",
    outcomes    = "binary requires exactly two levels (a 3-level factor is an error, not a 0/1/2 recode); bare proportions are refused because a proportion without its denominator has no defined binomial variance",
    fit_stat    = "'Deviance explained' = 1 - deviance/null.deviance. For binomial on 0/1 this equals McFadden's R^2 exactly (the saturated log-likelihood is zero); it does not for other families, which is why the display does not name a pseudo-R^2 family.",
    performance = "ROC, AUC and accuracy are APPARENT -- evaluated on the fitting sample, not cross-validated",
    caveat      = "not runtime-verified in the audit environment: refund would not install for R 4.3.3 there"
  ),

  run_functional_clustering = list(
    file        = "server/60_clustering.R",
    orientation = "data_matrix is SUBJECTS x TIME",
    estimand    = "a partition of subjects into k groups",
    geometry    = "cluster centres for every objective statistic come from data_matrix, the SAME matrix the assignments were computed on. Taking centres from values$fd_obj (always original scale) while members came from a standardised data_matrix made WCSS, R^2 and Calinski-Harabasz undefined and could drive between-SS negative.",
    display     = "cluster_means_raw holds the original-units means for plotting when standardisation is on",
    determinism = "seeded restarts (set.seed(100 + restart))",
    missing     = "complete cases"
  )
)

# Look a contract up from the console: fck_contract("fit_cosinor_nonlinear")
fck_contract <- function(name) {
  if (!name %in% names(FCK_CONTRACTS))
    stop("No contract recorded for '", name, "'. Have: ",
         paste(names(FCK_CONTRACTS), collapse = ", "))
  x <- FCK_CONTRACTS[[name]]
  cat("\n", name, "\n", strrep("-", nchar(name)), "\n", sep = "")
  for (f in names(x)) cat(sprintf("  %-12s %s\n", f, x[[f]]))
  invisible(x)
}
