# ==============================================================================
# server/08n_helpers_equations.R -- the equation behind each reported number
#
# Every quantity the app reports has an (i) beside it that opens the formula
# actually computed, written for the model that was fitted: the trend term is
# the one chosen, the harmonics the ones fitted, the log scale stated when it is
# on. One dictionary, so a formula is written once and every panel shows the
# same one. The formulas restate the code (dance_trend_factors and friends);
# they are not a second implementation.
# ==============================================================================

if (!exists("%||%", mode = "function")) `%||%` <- function(a, b) if (is.null(a)) b else a

# context: trend ("none", "linear", "log", "exp_sat"), K, period, approach
# ("mixed" / "two_stage"), log (TRUE when fitted on ln(y + c)), offset c,
# tau_mode ("fixed" / "estimate")
dance_eq_context <- function(mod = NULL) {
  if (is.null(mod)) return(list(trend = "none", K = 1L, period = 24, approach = "mixed", log = FALSE, offset = 0))
  sc <- mod$scale %||% list(type = "original", offset = 0)
  list(trend = mod$trend_type %||% "none", K = as.integer(mod$n_harmonics %||% 1L), period = mod$period %||% 24,
       approach = mod$approach %||% "mixed", log = identical(sc$type, "log"), offset = sc$offset %||% 0,
       tau_mode = mod$tau_spec$mode %||% "fixed", T = diff(range(mod$time_vec %||% c(0, 0))))
}

.eq <- function(...) paste0('<div style="font-family:Georgia,serif;font-size:15px;margin:8px 0 10px 0;',
                            'padding:8px 12px;background:#f6f8fa;border-left:3px solid #1f5f7a">', ..., '</div>')
.p <- function(...) paste0("<p>", ..., "</p>")
.b <- function(i, s = NULL) paste0("<i>&beta;</i><sub>", i, if (!is.null(s)) paste0(",", s), "</sub>")

dance_eq_trend_term <- function(trend) switch(trend,
  linear  = "<i>f</i>(<i>t</i>) = <i>b</i>&middot;<i>t</i>",
  log     = "<i>f</i>(<i>t</i>) = <i>b</i>&middot;ln(1 + <i>t</i>)",
  exp_sat = "<i>f</i>(<i>t</i>) = <i>A</i>&middot;(1 &minus; e<sup>&minus;<i>t</i>/<i>&tau;</i></sup>)",
  "<i>f</i>(<i>t</i>) = 0 (no trend)")

dance_eq_model <- function(ctx) {
  harm <- if (ctx$K >= 1L) paste0(" + &Sigma;<sub><i>h</i>=1</sub><sup>", ctx$K, "</sup> [", .b("c", "<i>h</i>"),
                                  " cos(2&pi;<i>h t</i>/<i>P</i>) + ", .b("s", "<i>h</i>"), " sin(2&pi;<i>h t</i>/<i>P</i>)]") else ""
  y <- if (isTRUE(ctx$log)) paste0("ln(<i>y</i>(<i>t</i>)", if (ctx$offset != 0) paste0(" + ", format(ctx$offset)), ")") else "<i>y</i>(<i>t</i>)"
  paste0(.eq(y, " = ", .b(0), " + <i>f</i>(<i>t</i>)", harm, " + <i>&epsilon;</i>"),
         .p(dance_eq_trend_term(ctx$trend), "; <i>t</i> in hours from the first sampling time, <i>P</i> = ",
            format(ctx$period), " h."),
         .p(if (identical(ctx$approach, "mixed"))
              "Mixed-effects: one model for all observations; each design cell has its own coefficients and each participant a random deviation from them."
            else "Two-stage: the model is fitted to each participant; the group values are the means of the participants' values."))
}

# key: an outcome key as the tables use it (fitted_mean, mesor, change,
# amplitude1, phase2, rhythm1, full, ...) or a test of one ("test:<key>").
dance_equation <- function(key, ctx = dance_eq_context()) {
  # zero:<key> -- a single trajectory's test against zero
  if (startsWith(key, "zero:")) {
    k0 <- sub("^zero:", "", key)
    what <- if (k0 == "circadian") paste0("every harmonic coefficient (", .b("c", "1"), ", ", .b("s", "1"), ", &hellip;)")
            else if (grepl("^rhythm[0-9]+$", k0)) { h0 <- sub("^rhythm", "", k0); paste0("(", .b("c", h0), ", ", .b("s", h0), ")") }
            else "the trend coefficient (the change over the window)"
    return(list(title = "Test against zero (one trajectory)", html = paste0(
      .eq("H<sub>0</sub>: ", what, " = 0"),
      .p(if (identical(ctx$approach, "mixed"))
           "A Wald F test of those fixed effects against zero (Satterthwaite or Kenward-Roger df): is there a rhythm, or a trend, at all?"
         else paste0("Population-mean cosinor (Bingham et al. 1982): a one-sample Hotelling <i>T</i><sup>2</sup> on the participants' coefficients, ",
                     "<i>F</i> = (<i>n</i> &minus; <i>p</i>)<i>T</i><sup>2</sup>/[<i>p</i>(<i>n</i> &minus; 1)] on <i>p</i> and <i>n</i> &minus; <i>p</i> df; ",
                     "for the trend, a one-sample <i>t</i> test of the participants' changes.")))))
  }
  test <- startsWith(key, "test:"); k <- sub("^test:", "", key)
  # a harmonic number only on the per-harmonic keys (rate0 and t50 end in digits too)
  harm_key <- grepl("^(amplitude|phase|rhythm)[0-9]+$", k)
  h <- if (harm_key) as.integer(sub("^[a-z]+([0-9]+)$", "\\1", k)) else NA_integer_
  base <- if (harm_key) sub("[0-9]+$", "", k) else k
  win <- if (is.finite(ctx$T %||% NA) && ctx$T > 0) sprintf("<i>T</i> = %s h, the observed window", format(round(ctx$T, 2))) else "<i>T</i> = the observed window"
  mean_of_f <- switch(ctx$trend,
    linear  = "<i>b</i>&middot;<i>T</i>/2",
    log     = "<i>b</i>&middot;[(<i>T</i> + 1) ln(1 + <i>T</i>) &minus; <i>T</i>]/<i>T</i>",
    exp_sat = "<i>A</i>&middot;[1 &minus; (<i>&tau;</i>/<i>T</i>)(1 &minus; e<sup>&minus;<i>T</i>/<i>&tau;</i></sup>)]",
    "0")
  chg <- switch(ctx$trend,
    linear = "<i>b</i>&middot;<i>T</i>", log = "<i>b</i>&middot;ln(1 + <i>T</i>)",
    exp_sat = "<i>A</i>&middot;(1 &minus; e<sup>&minus;<i>T</i>/<i>&tau;</i></sup>)", "0")
  rate0 <- switch(ctx$trend, linear = "<i>b</i>", log = "<i>b</i>", exp_sat = "<i>A</i>/<i>&tau;</i>", "0")
  logn <- if (isTRUE(ctx$log)) .p("Fitted on the log scale: read back on the original scale as e<sup><i>x</i></sup> &minus; <i>c</i>",
                                   if (base %in% c("fitted_mean", "mesor", "intercept")) " (a geometric mean when <i>c</i> = 0)"
                                   else if (base == "amplitude") "; an amplitude <i>A</i> as the peak-to-trough ratio e<sup>2<i>A</i></sup>"
                                   else if (base %in% c("change", "sat_change")) "; a change as the ratio e<sup>change</sup>" else "",
                                   ".") else ""
  grp <- if (identical(ctx$approach, "mixed"))
    .p("Per cell: a linear combination <i>L</i>&prime;<i>&beta;&#770;</i> of the cell's fixed effects, with interval <i>L</i>&prime;<i>&beta;&#770;</i> &plusmn; <i>t</i>&middot;SE (Satterthwaite or Kenward-Roger df).")
  else .p("Per group: the mean of the participants' values, with interval mean &plusmn; <i>t</i><sub><i>n</i>&minus;1</sub>&middot;SD/&radic;<i>n</i>.")
  out <- switch(if (test) paste0("test_", base) else base,
    full = , test_full = list("Complete fitted trajectory", paste0(dance_eq_model(ctx),
      .p("<b>Test</b>, H<sub>0</sub>: every coefficient (", .b(0), ", the trend, every ", .b("c", "<i>h</i>"), " and ",
         .b("s", "<i>h</i>"), ") is the same in every group. ",
         if (identical(ctx$approach, "mixed")) "A joint Wald F test of all the group contrasts at once."
         else "A one-way MANOVA (Wilks' &Lambda;) on the participants' coefficient vectors."))),
    fitted_mean = list("Time-averaged fitted response", paste0(
      .eq("<i>&#563;&#772;</i> = (1/<i>T</i>) &int;<sub>0</sub><sup><i>T</i></sup> <i>&#375;</i>(<i>t</i>) d<i>t</i>"),
      .p("The fitted curve, rhythm included, averaged over ", win, ". Over whole cycles the harmonics average to zero and this equals the time-averaged midline level."),
      grp, logn)),
    mesor = list("Time-averaged midline level", paste0(
      .eq("<i>M&#772;</i> = (1/<i>T</i>) &int;<sub>0</sub><sup><i>T</i></sup> [", .b(0), " + <i>f</i>(<i>t</i>)] d<i>t</i> = ",
          .b(0), " + ", mean_of_f),
      .p("The midline (intercept + trend, harmonics excluded) averaged over ", win,
         ". It generalizes the classical MESOR to a time-varying baseline; without a trend it is ", .b(0), "."),
      grp, logn)),
    intercept = list("Initial midline (t = 0)", paste0(
      .eq("<i>m</i>(0) = ", .b(0), " + <i>f</i>(0) = ", .b(0)),
      .p("The midline at the first sampling time, where every trend term is zero."), grp, logn)),
    change = , sat_change = list("Total trend change over the window", paste0(
      .eq("&Delta; = <i>f</i>(<i>T</i>) &minus; <i>f</i>(0) = ", chg), .p(dance_eq_trend_term(ctx$trend), "; ", win, "."),
      grp, logn)),
    avg_rate = , sat_avg_rate = list("Average trend rate", paste0(
      .eq("&Delta;/<i>T</i> = [<i>f</i>(<i>T</i>) &minus; <i>f</i>(0)]/<i>T</i> = ", chg, "/<i>T</i>"),
      .p("The total change per hour of the window: the same quantity as the change in other units, so it shares its test."), grp)),
    rate0 = , init_slope = list("Initial trend rate", paste0(
      .eq("<i>f</i>&prime;(0) = ", rate0), .p("The slope of the trend at the first sampling time."), grp)),
    A_sat = list("Saturating rise A", paste0(.eq("<i>f</i>(<i>t</i>) = <i>A</i>(1 &minus; e<sup>&minus;<i>t</i>/<i>&tau;</i></sup>) &rarr; <i>A</i> as <i>t</i> &rarr; &infin;"),
      .p("The asymptote of the trend above the intercept."), grp)),
    tau = list("Time constant tau", paste0(.eq("<i>f</i>(<i>&tau;</i>) = <i>A</i>(1 &minus; e<sup>&minus;1</sup>) &asymp; 0.632&middot;<i>A</i>"),
      .p(if (identical(ctx$tau_mode, "fixed")) "Held at the value set, not estimated." else "Estimated (profile likelihood or per participant, as the run states)."))),
    t50 = list("Half-rise time", .eq("<i>t</i><sub>50</sub> = <i>&tau;</i>&middot;ln 2")),
    r_squared = list("R-squared (per participant)", paste0(
      .eq("<i>R</i><sup>2</sup> = 1 &minus; SS<sub>res</sub> / SS<sub>tot</sub>"),
      .p("Of this participant's fit: residual against total sum of squares about their mean. Undefined when there is no variation."))),
    pct_rhythm = list("Percent of variance explained", paste0(.eq("100&middot;<i>R</i><sup>2</sup>"),
      .p("The whole fitted model (trend included), in percent; not the rhythm's share alone."))),
    zero_amp = list("p-value: no rhythm (per participant)", paste0(
      .eq("<i>F</i> = [(SS<sub>res,trend only</sub> &minus; SS<sub>res</sub>) / 2<i>K</i>] / [SS<sub>res</sub> / (<i>n</i> &minus; <i>p</i>)]"),
      .p("H<sub>0</sub>: every harmonic's amplitude is zero given the trend: the full fit against the fit with the trend alone, on 2<i>K</i> and <i>n</i> &minus; <i>p</i> df (<i>p</i> = parameters of the full model)."))),
    trend = , test_trend = list("Trend coefficient", paste0(.eq(dance_eq_trend_term(ctx$trend)),
      .p("The coefficient of the non-periodic trend (", switch(ctx$trend, linear = "<i>b</i>, per hour", log = "<i>b</i>",
         exp_sat = "<i>A</i>, the rise above the intercept", "none"), "); ",
         if (startsWith(key, "test:")) "tested as an F test of its cell contrasts." else "per cell from the fixed effects."))),
    trend_linear = list("Linear trend coefficient", .eq("<i>f</i>(<i>t</i>) = <i>b</i>&middot;<i>t</i>: <i>b</i> per hour")),
    trend_log = list("Logarithmic trend coefficient", .eq("<i>f</i>(<i>t</i>) = <i>b</i>&middot;ln(1 + <i>t</i>)")),
    amplitude = list(sprintf("Amplitude of harmonic %d", h), paste0(
      .eq("<i>A</i><sub>", h, "</sub> = &radic;(", .b("c", h), "<sup>2</sup> + ", .b("s", h), "<sup>2</sup>)"),
      .p("Half the peak-to-trough excursion of harmonic ", h, " (period ", format(ctx$period / h), " h) on its own. ",
         if (identical(ctx$approach, "mixed")) "From the cell's fixed effects." else
           "From the group's mean coefficient vector (population-mean cosinor), interval from Bingham et al.'s joint region."), logn)),
    phase = list(sprintf("Acrophase of harmonic %d", h), paste0(
      .eq("<i>&phi;</i><sub>", h, "</sub> = atan2(", .b("s", h), ", ", .b("c", h), "),&nbsp;&nbsp; clock time = <i>t</i><sub>0</sub> + <i>&phi;</i><sub>",
          h, "</sub>&middot;<i>P</i>/(2&pi;&middot;", h, ") (mod <i>P</i>/", h, ")"),
      .p("The time of the first maximum of harmonic ", h, " alone, not of the complete fitted curve. Not estimable when the joint (cos, sin) confidence region contains the origin."))),
    test_circadian = , circadian = list("Rhythmic component, all harmonics jointly", paste0(
      .eq("H<sub>0</sub>: (", .b("c", "1"), ", ", .b("s", "1"), ", &hellip;, ", .b("c", "<i>K</i>"), ", ", .b("s", "<i>K</i>"),
          ") equal in every group"),
      .p("A 2<i>K</i>-df joint test of every harmonic's cosine and sine coefficients (here <i>K</i> = ", ctx$K, "): ",
         if (identical(ctx$approach, "mixed")) "a Wald F test of the cell contrasts." else "a MANOVA on the participants' coefficients.",
         " The overall rhythm difference; neither an amplitude nor an acrophase test alone."))),
    test_rhythm = , rhythm = list(sprintf("%s h rhythmic component, joint test", format(ctx$period / h)), paste0(
      .eq("H<sub>0</sub>: (", .b("c", h), ", ", .b("s", h), ") equal in every group"),
      .p("A 2-df joint test of the cosine and sine coefficients of harmonic ", h, ": ",
         if (identical(ctx$approach, "mixed")) "a Wald F test of the cell contrasts." else "a MANOVA (Wilks' &Lambda;) on the participants' pairs.",
         " It tests amplitude and acrophase together, without the instability of testing an angle."))),
    NULL)
  if (is.null(out) && test)
    out <- list("Between-group test", paste0(
      .eq("H<sub>0</sub>: the quantity is equal in every group"),
      .p(if (identical(ctx$approach, "mixed")) "An F test of the linear contrast of the cells' fixed effects (Satterthwaite or Kenward-Roger df)."
         else "A one-way ANOVA on the participants' values (F with <i>g</i> &minus; 1 and <i>n</i> &minus; <i>g</i> df)."),
      .p("Unadjusted: a secondary comparison, read below the omnibus test.")))
  if (is.null(out)) out <- list(k, .p("No equation recorded for this quantity."))
  list(title = out[[1]], html = out[[2]])
}

# DLMO and functional ANOVA, which have no model context
dance_equation_fixed <- function(key) switch(key,
  dlmo_fixed = list(title = "DLMO, fixed threshold", html = paste0(
    .eq("<i>t</i>* = <i>t</i><sub>1</sub> + (<i>&theta;</i> &minus; <i>y</i><sub>1</sub>)&middot;(<i>t</i><sub>2</sub> &minus; <i>t</i><sub>1</sub>)/(<i>y</i><sub>2</sub> &minus; <i>y</i><sub>1</sub>)"),
    .p("Linear interpolation between the two observed samples bracketing the first rise above the threshold <i>&theta;</i> (<i>y</i><sub>1</sub> &le; <i>&theta;</i> &lt; <i>y</i><sub>2</sub>), the following sample(s) also above it; only when the two samples are no further apart than the interpolation limit."))),
  dlmo_variable = list(title = "DLMO, variable threshold", html = paste0(
    .eq("<i>&theta;</i> = <i>&#563;&#772;</i><sub>1..<i>n</i></sub> + <i>k</i>&middot;SD(<i>y</i><sub>1..<i>n</i></sub>)"),
    .p("The mean plus <i>k</i> SD of the first <i>n</i> observed samples of the profile (Voultsios et al. 1997; default <i>n</i> = 3, <i>k</i> = 2); the onset then as for the fixed threshold."))),
  dlmo_hockey = list(title = "DLMO, hockey-stick", html = paste0(
    .p("Danilenko et al. (2014), computed by the dlmoR package: a horizontal baseline and an ascending segment fitted to the profile, the onset where they meet, searched on a coarse then a fine grid."))),
  phase_angle = list(title = "Phase angle of entrainment", html = paste0(
    .eq("<i>&psi;</i> = bedtime &minus; DLMO"), .p("In hours; positive when the onset precedes bedtime. The bedtime used is the one within 12 h of the onset."))),
  fanova_F = list(title = "Pointwise F statistic", html = paste0(
    .eq("<i>F</i>(<i>t</i>) = [SSB(<i>t</i>)/(<i>g</i> &minus; 1)] / [SSW(<i>t</i>)/(<i>N</i> &minus; <i>g</i>)]"),
    .p("At each of 100 evaluation points: between-group over within-group variation of the curves. Its p-value is by permutation of the group labels, then FDR-adjusted across the points."))),
  fanova_L2 = list(title = "Global L2 statistic", html = paste0(
    .eq("<i>L</i><sup>2</sup> = ( &int; [SSB(<i>t</i>)/<i>N</i>]<sup>2</sup> d<i>t</i> )<sup>1/2</sup>"),
    .p("One number for the whole curve; its p-value by permutation of the group labels."))),
  fanova_p = list(title = "Permutation p-value and FDR", html = paste0(
    .eq("<i>p</i> = (1 + #{<i>T</i>* &ge; <i>T</i>}) / (1 + <i>B</i>)"),
    .p("<i>B</i> permutations of the group labels; never 0, and not computed when no permutation gives a defined statistic. Pointwise p-values are then adjusted across the evaluation points by Benjamini-Hochberg: <i>p</i><sub>(<i>i</i>)</sub>&middot;<i>m</i>/<i>i</i>, made monotone."))),
  fanova_eta = list(title = "Pointwise eta-squared", html = .eq("<i>&eta;</i><sup>2</sup>(<i>t</i>) = SSB(<i>t</i>) / SST(<i>t</i>)")),
  NULL)

# The (i) link. Clicking it opens the equation; it does not also click the row
# it sits in.
dance_eq_link <- function(key, label = NULL) {
  shiny::tags$a(href = "#", title = "Show the equation", style = "margin-left:5px;color:#1f5f7a;text-decoration:none",
    onclick = sprintf(paste0("event.preventDefault();event.stopPropagation();",
                             "Shiny.setInputValue('dance_eq', {key: '%s', nonce: Math.random()}, {priority: 'event'});"), key),
    shiny::icon("info-circle"), if (!is.null(label)) label)
}
