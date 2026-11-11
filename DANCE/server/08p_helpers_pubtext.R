# ==============================================================================
# server/08p_helpers_pubtext.R -- publication-ready text for the cosinor
# ==============================================================================
# A Statistical Analysis section, a Results section and a Table 1, written in
# APA style from the SAME tables the app shows: the outcome table
# (dance_outcomes_mixed / dance_outcomes_two_stage: estimates and tests by
# question) and the pairwise table (every quantity's pairs, with the
# adjustment chosen on the comparison tab). Nothing is recomputed here and no
# number is typed in: each sentence is chosen by a rule from the tests
# (significant at alpha or not) and the group order from the estimates.
#
#   o     the outcome table (list with estimates, tests, period, clock_origin,
#         approach, trend, structure)
#   pw    the pairwise table (harmonic_pairwise_all(): quantity, comparison
#         "a vs b" = a minus b, estimate, ci_lo, ci_hi, p, p_adjusted, adjustment)
#   info  the run's settings: approach, factor, roles, n_part, n_obs, n_time,
#         n_missing, window, period, K, trend, tau, tau_mode, re_terms,
#         re_uncorrelated, residuals, phi, df_method, adjust, r2_marginal,
#         r2_conditional, r2_ts (mean, min, max), group_n (named), versions,
#         origin_text, log_scale, units
# Returns markdown lines.
dance_pub_cosinor <- function(o, pw, info, alpha = 0.05) {
  `%|%` <- function(a, b) if (is.null(a) || length(a) == 0L || (length(a) == 1L && is.na(a))) b else a
  L <- character(0); add <- function(...) L <<- c(L, paste0(...)); blank <- function() L <<- c(L, "")
  es <- o$estimates; ts <- o$tests
  mixed <- identical(info$approach, "mixed")
  f2 <- function(x) if (!is.finite(x)) "NA" else formatC(x, format = "f", digits = 2)
  # "a and b"; "a, b, and c"
  and_list <- function(x, conj = "and") if (length(x) <= 1L) x else if (length(x) == 2L) paste(x, collapse = paste0(" ", conj, " ")) else
    paste0(paste(x[-length(x)], collapse = ", "), ", ", conj, " ", x[length(x)])
  pfmt <- function(p) if (!is.finite(p)) "p = NA" else if (p < .001) "*p* < .001" else
    paste0("*p* = ", sub("^0", "", formatC(p, format = "f", digits = 3)))
  Ftxt <- function(r) sprintf("*F*(%s, %s) = %s, %s", format(r$df1),
                              if (is.finite(r$df2)) formatC(r$df2, format = "f", digits = 1) else "inf",
                              f2(r$statistic), pfmt(r$p))
  cells <- unique(es$cell)
  multi <- length(cells) > 1L
  numeric_lv <- all(grepl("^[0-9.]+$", cells))
  gname <- function(cl) if (!multi) "the sample" else if (numeric_lv) paste("Group", cl) else sprintf("the %s group", cl)
  gname_pl <- function(a, b) if (numeric_lv) sprintf("Groups %s and %s", a, b) else sprintf("the %s and %s groups", a, b)
  cap <- function(x) paste0(toupper(substr(x, 1, 1)), substr(x, 2, nchar(x)))
  fac <- gsub("[._]+", " ", info$factor %|% "group")
  # "between age groups", "between the groups" (a factor called group)
  fac_groups <- if (identical(tolower(fac), "group")) "the groups"
                else if (grepl("groups$", fac, ignore.case = TRUE)) paste("the", fac)
                else if (grepl("group$", fac, ignore.case = TRUE)) paste0("the ", fac, "s")
                else paste("the", fac, "groups")
  est <- function(q, cl) { r <- es[es$outcome == q & es$cell == cl, , drop = FALSE]; if (nrow(r)) r[1, ] else NULL }
  val_ci <- function(q, cl) {
    r <- est(q, cl); if (is.null(r)) return("NA")
    if (grepl("^phase", q)) return(dance_outcome_value_text(r, o$period, o$clock_origin, digits = 2))
    if (is.finite(r$lo) && is.finite(r$hi)) sprintf("%s, 95%% CI [%s, %s]", f2(r$estimate), f2(r$lo), f2(r$hi))
    else f2(r$estimate)
  }
  tst <- function(q, dom = NULL) {
    # the first outcome named that the table has (q in order of preference)
    for (qq in q) {
      r <- ts[ts$outcome == qq & (is.null(dom) | ts$domain %in% dom), , drop = FALSE]
      if (nrow(r)) return(r[1, ])
    }
    NULL
  }
  sig <- function(r) !is.null(r) && is.finite(r$p) && r$p < alpha
  adj <- (pw$adjustment %|% "none")[1]
  padj_word <- switch(adj, holm = "Holm-adjusted", bonferroni = "Bonferroni-adjusted", "")
  # one quantity's pairwise sentences: significant pairs with the difference
  # oriented so that it is positive, the others listed with their p
  pair_txt <- function(q, noun, higher = "higher", lower = "lower") {
    if (is.null(pw)) return("")
    d <- pw[pw$quantity == q & is.finite(pw$estimate), , drop = FALSE]
    if (!nrow(d)) return("")
    pp <- if (identical(adj, "none")) d$p else d$p_adjusted
    parts <- strsplit(d$comparison, " vs ", fixed = TRUE)
    s_txt <- character(0); ns_txt <- character(0)
    for (i in seq_len(nrow(d))) {
      a <- parts[[i]][1]; b <- parts[[i]][2]
      e <- d$estimate[i]; lo <- d$ci_lo[i]; hi <- d$ci_hi[i]
      if (e < 0) { tmp <- a; a <- b; b <- tmp; e <- -e; l2 <- -hi; hi <- -lo; lo <- l2 }
      if (is.finite(pp[i]) && pp[i] < alpha)
        s_txt <- c(s_txt, sprintf("%s had a %s %s than %s (difference = %s%s, %s)", gname(a), higher, noun,
                                  gname(b), f2(e), if (is.finite(lo)) sprintf(", 95%% CI [%s, %s]", f2(lo), f2(hi)) else "",
                                  pfmt(pp[i])))
      else ns_txt <- c(ns_txt, sprintf("%s (%s)", gname_pl(parts[[i]][1], parts[[i]][2]), pfmt(pp[i])))
    }
    paste0(if (length(s_txt)) paste0("Pairwise comparisons", if (nzchar(padj_word)) paste0(" (", padj_word, ")") else "",
                                     " indicated that ", paste(s_txt, collapse = "; "), ". ") else "",
           if (length(ns_txt)) paste0(if (length(s_txt)) "No significant difference was observed between " else
                                        paste0("Pairwise comparisons", if (nzchar(padj_word)) paste0(" (", padj_word, ")") else "",
                                               " showed no significant difference between "),
                                      and_list(ns_txt, "or"), ".") else "")
  }
  ordered_groups <- function(q, decreasing = TRUE) {
    v <- vapply(cells, function(cl) { r <- est(q, cl); if (is.null(r)) NA_real_ else r$estimate }, numeric(1))
    cells[order(v, decreasing = decreasing)]
  }
  trend <- o$trend %|% "none"; K <- as.integer(info$K %|% 0L)
  has_trend <- !identical(trend, "none")
  trend_word <- switch(trend, linear = "linear", log = "logarithmic", exp_sat = "saturating exponential", "")
  q_level <- if (any(es$outcome == "fitted_mean")) "fitted_mean" else "mesor"
  q_change <- if (any(es$outcome == "sat_change")) "sat_change" else "change"
  win <- info$window %|% "the observation window"

  # ---------------------------------------------------------------- METHODS
  add("### Statistical Analysis"); blank()
  gn <- info$group_n
  add("Longitudinal trajectories were analysed using ", if (mixed) "mixed-effects" else "two-stage", " harmonic (cosinor) regression. ",
      "The dataset comprised ", info$n_part, " participants, each providing up to ", info$n_time,
      " observations over ", win, ", yielding ", format(info$n_obs, big.mark = ","), " observations",
      if (isTRUE(info$n_missing == 0)) " without missing values" else sprintf(" (%s missing)", format(info$n_missing, big.mark = ",")), ".",
      if (multi && length(gn)) paste0(" Participants were classified by ", fac, " (",
                                      paste(sprintf("%s: *n* = %d", if (numeric_lv) paste("Group", names(gn)) else names(gn), as.integer(gn)), collapse = "; "), ").") else "")
  blank()
  harm_txt <- if (K >= 1L) paste0(if (K == 1L) sprintf("a %s-hour harmonic component", format(info$period))
                                  else sprintf("%d harmonics of a %s-hour period (%s h)", K, format(info$period),
                                               paste(format(round(info$period / seq_len(K), 2)), collapse = ", ")),
                                  ", represented by sine and cosine terms,") else NULL
  trend_txt <- if (has_trend) paste0("a non-periodic ", trend_word, " trend",
    if (identical(trend, "exp_sat")) {
      if (identical(info$tau_mode, "fixed") && !mixed) sprintf(" with its time constant held at \u03c4 = %s h for every participant", format(signif(info$tau, 3)))
      else if (identical(info$tau_mode, "estimate") && !mixed) " with its time constant estimated for each participant"
      else if (identical(info$tau_mode, "fixed")) sprintf(" with a fixed time constant (τ = %s h), so that the magnitude of the increase could vary between groups while its temporal shape was common", format(signif(info$tau, 3)))
      else if (identical(info$tau_mode, "per_cell")) " with its magnitude and time constant estimated in each group (nonlinear mixed model)"
      else sprintf(" with a time constant estimated by profile likelihood (τ = %s h)", format(signif(info$tau, 3)))
    } else "") else NULL
  comps <- c(harm_txt, trend_txt)
  if (mixed) {
    add("A linear mixed-effects model was fitted to the raw observations using restricted maximum likelihood (REML). ",
        "The model included ", paste(comps, collapse = " and "), ".",
        if (!is.null(info$origin_text)) paste0(" ", info$origin_text) else "",
        if (multi) paste0(" ", cap(fac), " was included as a categorical ", info$role %|% "between-subjects",
                          " factor, interacting with the intercept", if (has_trend) ", trend" else "",
                          if (K >= 1L) " and harmonic terms" else "", ".") else "",
        " Participant-specific random effects were specified for the ",
        if (length(info$re_terms) > 1L) paste0(paste(info$re_terms[-length(info$re_terms)], collapse = ", "), " and ",
                                               info$re_terms[length(info$re_terms)]) else info$re_terms,
        if (isTRUE(info$re_uncorrelated)) ", assumed uncorrelated." else ", with an unstructured covariance matrix.",
        if (identical(info$residuals, "car1"))
          sprintf(" Residuals within a curve were modelled as a continuous-time first-order autoregressive process (φ = %s per hour).", formatC(info$phi, format = "f", digits = 2))
        else " Residuals were assumed to be independent and homoscedastic.",
        if (isTRUE(info$log_scale)) " The response was analysed on the natural-log scale." else "")
  } else {
    add("A cosinor model was fitted to each participant's series by ", if (identical(trend, "exp_sat")) "nonlinear " else "",
        "least squares. It comprised ", paste(comps, collapse = " and "), ".",
        if (!is.null(info$origin_text)) paste0(" ", info$origin_text) else "",
        " The participants' estimates were then summarised by group: group amplitudes and acrophases by the population-mean cosinor (Bingham et al., 1982).",
        if (isTRUE(info$log_scale)) " The response was analysed on the natural-log scale." else "")
  }
  blank()
  add("The model provided estimates of the time-averaged response level",
      if (has_trend) ", the non-periodic change across the observation period" else "",
      if (K >= 1L) ", and the amplitude and acrophase of each harmonic component" else "", ". ",
      if (multi) paste0(
        if (mixed) paste0("Group differences were evaluated using omnibus and component-specific *F* tests with ",
                          info$df_method, "-adjusted denominator degrees of freedom. Pairwise group comparisons were conducted using model-based contrasts. ")
        else paste0("Group differences were evaluated with a multivariate test of the participants' coefficient vectors (MANOVA, Wilks' Λ), ",
                    "one-way ANOVAs of the per-participant estimates, and the population-mean cosinor for amplitude and acrophase. ",
                    "Pairwise comparisons used Welch's *t* tests, and Watson-Williams tests for acrophases. "),
        if (identical(adj, "none")) "These comparisons were not adjusted for multiple testing. "
        else sprintf("*p* values were %s within each set of pairwise comparisons. ", tolower(padj_word))) else
        "Each component was tested against zero. ",
      "Parameter estimates are presented with 95% confidence intervals (CIs)",
      if (mixed) ", obtained using Wald-based methods for linear contrasts and simulation-based methods for harmonic amplitude and acrophase." else ".")
  blank()
  add("Analyses were conducted in R (version ", info$versions[["R"]] %|% "", "), using the DANCE analysis environment",
      if (mixed) paste0(" and the *lme4* (", info$versions[["lme4"]] %|% "", ") and *lmerTest* (", info$versions[["lmerTest"]] %|% "", ") packages.")
      else ".")
  blank()

  # ---------------------------------------------------------------- RESULTS
  add("### Results"); blank()
  add("#### Overall model fit", if (multi) paste0(" and ", fac, " differences") else ""); blank()
  fit_txt <- if (mixed && is.finite(info$r2_marginal %|% NA))
    sprintf("The mixed-effects harmonic regression model had a marginal *R*² of %s and a conditional *R*² of %s: the fixed effects explained %s%% of the variance, and %s%% when participant-specific random effects were included.",
            sub("^0", "", formatC(info$r2_marginal, format = "f", digits = 3)), sub("^0", "", formatC(info$r2_conditional, format = "f", digits = 3)),
            formatC(100 * info$r2_marginal, format = "f", digits = 1), formatC(100 * info$r2_conditional, format = "f", digits = 1))
  else if (!mixed && length(info$r2_ts) == 3L)
    sprintf("The participants' cosinor fits had a mean *R*² of %s (range %s to %s).",
            sub("^0", "", formatC(info$r2_ts[1], format = "f", digits = 3)), sub("^0", "", formatC(info$r2_ts[2], format = "f", digits = 3)),
            sub("^0", "", formatC(info$r2_ts[3], format = "f", digits = 3)))
  else ""
  om <- tst("full", "omnibus")
  dom_tests <- list(level = tst(q_level, "level"), change = if (has_trend) tst(q_change, "accumulation"),
                    rhythm = if (K >= 1L) tst(c("circadian", "rhythm1"), "rhythm"))
  dom_word <- c(level = "overall response level", change = "the non-periodic change", rhythm = "the rhythmic component")
  sg <- names(dom_tests)[vapply(dom_tests, sig, logical(1))]
  nsg <- setdiff(names(dom_tests)[!vapply(dom_tests, is.null, logical(1))], sg)
  add(fit_txt,
      if (multi && !is.null(om)) paste0(
        " The omnibus test ", if (sig(om)) "revealed a significant" else "did not reveal a significant", " effect of ", fac,
        " on the overall fitted trajectories, ", Ftxt(om), ".",
        if (length(sg)) paste0(" Component-specific tests indicated that ", if (sig(om)) "these " else "", "differences were associated with ",
                               paste(dom_word[sg], collapse = " and "),
                               if (length(nsg)) paste0(", rather than ", paste(dom_word[nsg], collapse = " or ")) else "", ".")
        else " None of the component-specific tests was significant.") else "")
  blank()

  # overall response level
  if (multi) {
    add("#### Overall response level"); blank()
    r <- dom_tests$level; og <- ordered_groups(q_level)
    lab_lv <- if (q_level == "fitted_mean") "time-averaged fitted response" else "time-averaged midline level"
    add("The ", lab_lv, if (sig(r)) " differed significantly" else " did not differ significantly", " between ", fac_groups,
        if (!is.null(r)) paste0(", ", Ftxt(r)) else "", ". ",
        cap(sprintf("%s exhibited the highest estimated level (%s), followed by ", gname(og[1]), val_ci(q_level, og[1]))),
        paste(vapply(og[-1], function(cl) sprintf("%s (%s)", gname(cl), val_ci(q_level, cl)), ""), collapse = " and "), ".")
    blank()
    add(pair_txt(q_level, "time-averaged level")); blank()
    if (has_trend) {
      ri <- tst("intercept", "level")
      if (!is.null(ri)) {
        add(cap(fac_groups), " ", if (sig(ri)) "also differed" else "did not differ significantly",
            " in their estimated initial midline (the midline at ", info$origin_clock %|% "the first observation", "), ", Ftxt(ri), ". ",
            pair_txt("intercept", "initial midline"))
        blank()
      }
    }
  }
  # non-periodic change
  if (has_trend) {
    add("#### Non-periodic temporal change"); blank()
    r <- dom_tests$change; og <- ordered_groups(q_change)
    if (multi) {
      add(if (sig(r)) "Significant" else "No significant", " differences between ", fac_groups, " were observed in the non-periodic trend component",
          if (!is.null(r)) paste0(", ", Ftxt(r)) else "", ". The estimated change over ", win, " was greatest in ",
          gname(og[1]), " (", val_ci(q_change, og[1]), "), followed by ",
          paste(vapply(og[-1], function(cl) sprintf("%s (%s)", gname(cl), val_ci(q_change, cl)), ""), collapse = " and "), ".")
      blank()
      add(pair_txt(q_change, "change")); blank()
    } else {
      add("The estimated non-periodic change over ", win, " was ", val_ci(q_change, cells[1]),
          if (!is.null(r)) paste0(" (test against zero: ", Ftxt(r), ")") else "", ".")
      blank()
    }
    if (any(es$outcome == "change_origin")) {
      add("Measured from the time origin (", info$origin_clock %|% "", "), the change was ",
          paste(vapply(cells, function(cl) sprintf("%s in %s", val_ci("change_origin", cl), gname(cl)), ""), collapse = ", "),
          "; this extends the fitted trend before the first observation and is an extrapolation.")
      blank()
    }
    if (identical(trend, "exp_sat") && identical(info$tau_mode, "fixed"))
      { add("As the exponential time constant was fixed at ", format(signif(info$tau, 3)), " hours across groups, these differences reflect variation in the magnitude of the non-periodic increase rather than in its temporal rate constant."); blank() }
  }
  # rhythm
  if (K >= 1L) {
    add("#### Rhythmic component"); blank()
    r <- dom_tests$rhythm
    # a conclusion only where there is a test to state it with
    if (multi && !is.null(r)) add("The joint test of the harmonic coefficients ", if (sig(r)) "revealed" else "did not reveal",
                   " a significant overall difference in the rhythmic component between ", fac_groups, ", ", Ftxt(r), ".")
    else if (!is.null(r)) add("The rhythmic component was ", if (sig(r)) "" else "not ", "significantly different from zero, ", Ftxt(r), ".")
    blank()
    for (h in seq_len(K)) {
      qa <- sprintf("amplitude%d", h); qp <- sprintf("phase%d", h)
      av <- vapply(cells, function(cl) { r <- est(qa, cl); if (is.null(r)) "NA" else
        sprintf("%s%s in %s", f2(r$estimate), if (is.finite(r$lo)) sprintf(" (95%% CI [%s, %s])", f2(r$lo), f2(r$hi)) else "", gname(cl)) }, "")
      add("The estimated ", if (K > 1L) sprintf("H%d (%s h) ", h, format(round(info$period / h, 2))) else "", "amplitude was ",
          and_list(av),
          ". Corresponding acrophases occurred at ",
          { pv <- vapply(cells, function(cl) sprintf("%s in %s",
              sub("^(\\S+) \\[(.*)\\]$", "\\1 (95% CI [\\2])", val_ci(qp, cl)), gname(cl)), "")
            and_list(pv) },
          ".")
      blank()
      if (multi && !is.null(pw)) {
        da <- pw[pw$quantity == qa & is.finite(pw$p), , drop = FALSE]
        pa <- if (identical(adj, "none")) da$p else da$p_adjusted
        dp <- pw[pw$quantity == qp & is.finite(pw$p), , drop = FALSE]
        pph <- if (identical(adj, "none")) dp$p else dp$p_adjusted
        none_txt <- function(pp) sub("= ", "\u2265 ", pfmt(min(pp)))
        if (nrow(da) && nrow(dp) && !any(pa < alpha) && !any(pph < alpha)) {
          add(sprintf("No significant pairwise differences were observed in amplitude (all %s) or acrophase (all %s).",
                      none_txt(pa), none_txt(pph)))
          blank()
          next
        }
        if (nrow(da) && !any(pa < alpha)) {
          add(sprintf("No significant pairwise differences were observed in amplitude (all %s).", none_txt(pa)))
          blank()
        }
        if (nrow(da) && any(pa < alpha)) {
          add("Exploratory pairwise comparisons indicated ", pair_txt(qa, "amplitude", "greater", "smaller") |>
                sub(pattern = "^Pairwise comparisons( \\([^)]*\\))? indicated that ", replacement = ""),
              if (!sig(r) || identical(adj, "none")) paste0(" Given ", if (!sig(r)) "the nonsignificant omnibus test" else "",
                                                           if (!sig(r) && identical(adj, "none")) " and " else "",
                                                           if (identical(adj, "none")) "the absence of multiple-comparison adjustment" else "",
                                                           ", this finding should be interpreted cautiously.") else "")
          blank()
        }
        if (nrow(dp)) {
          add(if (any(pph < alpha)) pair_txt(qp, "acrophase (difference in hours)", "later", "earlier")
              else sprintf("No significant pairwise differences were observed in acrophase (all %s).",
                           none_txt(pph)))
          blank()
        }
      }
    }
  }
  # summary
  if (multi) {
    add("#### Summary of findings"); blank()
    hi_lv <- ordered_groups(q_level)[1]
    add("Taken together, ", if (length(sg)) paste0(fac, " differences in the fitted trajectories were predominantly associated with ",
                                                   paste(dom_word[sg], collapse = " and "), ".")
        else paste0("the analyses provided no significant evidence of ", fac, " differences in the components of the fitted trajectories."),
        if ("level" %in% sg) paste0(" ", cap(sprintf("%s exhibited the highest average response level.", gname(hi_lv)))) else "",
        if ("change" %in% sg) paste0(" ", cap(sprintf("%s showed the largest non-periodic change.", gname(ordered_groups(q_change)[1])))) else "",
        if (K >= 1L && length(sg) && !"rhythm" %in% sg) paste0(" There was no statistically significant overall evidence of ", fac, " differences in the rhythmic component.") else "")
    blank()
  }
  add("These findings describe variation in the periodic and non-periodic components of the observed trajectories and should not be interpreted as direct evidence of differences in underlying endogenous circadian or homeostatic processes.")
  blank()

  # ---------------------------------------------------------------- TABLE 1
  add("#### Table 1"); blank()
  add("*", if (mixed) "Model-based parameter estimates" else "Parameter estimates", if (multi) paste0(" by ", fac) else "", ".*"); blank()
  hdr <- vapply(cells, function(cl) sprintf("%s%s", if (numeric_lv) paste("Group", cl) else cl, if (!is.null(gn[[cl]])) sprintf(" (*n* = %d)", as.integer(gn[[cl]])) else ""), "")
  row <- function(lbl, q) paste0("| ", lbl, " | ", paste(vapply(cells, function(cl) {
    r <- est(q, cl); if (is.null(r)) return("")
    if (grepl("^phase", q)) return(dance_outcome_value_text(r, o$period, o$clock_origin, digits = 2))
    if (is.finite(r$lo)) sprintf("%s [%s, %s]", f2(r$estimate), f2(r$lo), f2(r$hi)) else f2(r$estimate)
  }, ""), collapse = " | "), " |")
  add("| Parameter | ", paste(hdr, collapse = " | "), " |")
  add("|---|", paste(rep("---|", length(cells)), collapse = ""))
  add(row(if (q_level == "fitted_mean") "Time-averaged response" else "Time-averaged midline level", q_level))
  if (has_trend) {
    add(row("Initial midline", "intercept"))
    add(row(sprintf("Non-periodic change (%s)", info$window_short %|% "over the window"), q_change))
    if (any(es$outcome == "change_origin")) add(row(sprintf("Change since the time origin (%s), extrapolated", info$origin_clock %|% ""), "change_origin"))
  }
  for (h in seq_len(K)) {
    lbl <- if (K > 1L) sprintf(" (H%d, %s h)", h, format(round(info$period / h, 2))) else sprintf(" (%s h)", format(info$period))
    add(row(paste0("Harmonic amplitude", lbl), sprintf("amplitude%d", h)))
    add(row(paste0("Harmonic acrophase", lbl), sprintf("phase%d", h)))
  }
  blank()
  same_lv <- all(vapply(cells, function(cl) {
    a <- est("fitted_mean", cl); b <- est("mesor", cl)
    !is.null(a) && !is.null(b) && abs(a$estimate - b$estimate) < 1e-6 }, logical(1)))
  add("*Note.* Values are ", if (mixed) "model-based estimates" else "group estimates from the per-participant fits",
      " with 95% confidence intervals.",
      if (has_trend) paste0(" The initial midline is the non-periodic component at ", info$origin_clock %|% "the first observation", ".") else "",
      if (has_trend) paste0(" Non-periodic change is the estimated change over ", win,
                            if (identical(trend, "exp_sat") && identical(info$tau_mode, "fixed"))
                              sprintf(", conditional on a fixed exponential time constant of %s hours", format(signif(info$tau, 3))) else "", ".") else "",
      if (K >= 1L) " Acrophase is the clock time of the maximum of the isolated harmonic component, not of the complete fitted trajectory." else "",
      if (K >= 1L && same_lv) " The time-averaged fitted response equals the time-averaged midline level because the harmonic component averages to zero over the complete observation period." else "",
      if (!mixed) " Group amplitudes and acrophases are population-mean cosinor estimates (Bingham et al., 1982)." else "")
  blank()
  sub("[ \t]+$", "", L)
}
