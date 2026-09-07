# ==============================================================================
# server/93_apa_report.R — the APA-style analysis report
#
# WHAT THIS IS FOR. Every other export writes numbers: a CSV of scores, a PDF of
# plots, an R script that recomputes them. None of them writes the SENTENCES a
# paper needs, and the gap between "the app printed F = 4.21" and "the Methods
# section says which estimator produced it, on what basis, with which
# multiplicity correction, at what permutation resolution" is where reporting
# errors live. This module writes both halves: the analytic choices, and the
# results, in APA 7th-edition style.
#
# THREE RULES IT FOLLOWS.
#
# 1. It reports only what was actually run. Every section is conditional on the
#    corresponding result existing in `values`. A report describing an analysis
#    you did not perform is worse than no report.
#
# 2. It carries the caveats with the numbers. This audit repeatedly found the
#    app claiming more than it had computed -- an "AIC" with no likelihood, a
#    "variance decomposition" that did not decompose, a standard error of zero
#    meaning "unavailable". Those corrections are worth nothing if the document
#    generated FOR PUBLICATION reintroduces the overstatement. So each section
#    ends with what its numbers do NOT establish, in the same voice.
#
# 3. It states the resolution of what it reports. A permutation p-value has a
#    floor at 1/(B+1). Reporting "p = .0099" without "B = 100, so p cannot go
#    below .0099" reports a bound as if it were an estimate.
#
# APA 7 formatting is implemented in fck_apa_* below: no leading zero on a
# quantity that cannot exceed 1, exact p to three decimals with "< .001" as the
# floor, statistics italicised, exact degrees of freedom.
# ==============================================================================


# ---- APA number formatting ---------------------------------------------------

# APA 7 sec. 6.36: no leading zero for a statistic that cannot exceed 1
# (p, r, R-squared, eta-squared); a leading zero for everything else.
fck_apa_num <- function(x, digits = 2, bounded = FALSE) {
  if (is.null(x) || length(x) == 0) return("--")
  x <- suppressWarnings(as.numeric(x))
  out <- vapply(x, function(z) {
    if (!is.finite(z)) return("--")
    s <- formatC(z, format = "f", digits = digits)
    if (bounded) s <- sub("^(-?)0\\.", "\\1.", s)
    s
  }, character(1))
  out
}

# APA 7 sec. 6.42: exact p to three decimals; "< .001" below that.
fck_apa_p <- function(p, digits = 3) {
  v <- fck_apa_pval(p, digits)
  if (identical(v, "--")) return("--")
  if (substr(v, 1, 1) == "<" || substr(v, 1, 1) == ">") paste0("*p* ", v)
  else paste0("*p* = ", v)
}

fck_apa_pval <- function(p, digits = 3) {
  if (is.null(p) || length(p) == 0) return("--")
  p <- suppressWarnings(as.numeric(p))[1]
  if (!is.finite(p)) return("--")
  if (p < .001) return("< .001")
  if (p > .999) return("> .999")
  sub("^0\\.", ".", formatC(p, format = "f", digits = digits))
}

fck_apa_F <- function(Fv, df1, df2, p = NULL) {
  if (!is.finite(suppressWarnings(as.numeric(Fv)))) return("*F* undefined")
  out <- sprintf("*F*(%s, %s) = %s", format(df1), format(df2), fck_apa_num(Fv, 2))
  if (!is.null(p)) out <- paste0(out, ", ", fck_apa_p(p))
  out
}

# AUDIT (P12.2): `bounded` was missing here, so the cosinor section printed
# "R-squared M = 0.958" -- a leading zero on a quantity that cannot exceed 1,
# which is an APA 7 error (6.36) in exactly the same way the p-values and ranges
# were. The SD of a bounded quantity is itself bounded, so both parts follow the
# rule. The APA test suite checked p-values and ranges and not this, which is
# why it shipped.
fck_apa_M <- function(x, digits = 2, bounded = FALSE) {
  x <- x[is.finite(x)]
  if (!length(x)) return("--")
  sprintf("*M* = %s, *SD* = %s", fck_apa_num(mean(x), digits, bounded),
          fck_apa_num(stats::sd(x), digits, bounded))
}

# AUDIT (P12.2): the pairwise controls carry internal ids (amplitude_1,
# acrophase_time_1, mesor_adj). Those were being printed into report prose --
# "Group differences in amplitude_1" -- which is not something that can go into
# a paper. One mapping, used by every sentence that names a parameter.
fck_cosinor_param_label <- function(p) {
  if (is.null(p) || !nzchar(p)) return("the selected parameter")
  if (grepl("^amplitude_([0-9]+)$", p))
    return(sprintf("the amplitude of harmonic %s", sub("^amplitude_", "", p)))
  if (grepl("^acrophase_time_([0-9]+)$", p))
    return(sprintf("the acrophase of harmonic %s", sub("^acrophase_time_", "", p)))
  switch(p,
    mesor          = "the constant term",
    mesor_adj      = "the MESOR (rhythm-adjusted mean)",
    value_at_start = "the predicted value at the first observation",
    p)
}

fck_cosinor_origin_label <- function(x) {
  if (is.null(x) || !nzchar(x)) return(NULL)
  switch(x,
    first_observation = "the first observation in each series",
    midnight          = "midnight",
    x)
}

# `bounded` propagates the no-leading-zero rule to both ends of a range: an
# eta-squared range printed as "0.007 to 0.241" is an APA error in the same way
# a single "0.007" would be.
fck_apa_range <- function(x, digits = 2, bounded = FALSE) {
  x <- x[is.finite(x)]
  if (!length(x)) return("--")
  sprintf("%s to %s", fck_apa_num(min(x), digits, bounded),
          fck_apa_num(max(x), digits, bounded))
}

# APA writes an eponymous procedure with its author's name capitalised. The app
# stores these as p.adjust() method strings ("bonferroni", "holm", "BH"), which
# are right as arguments and wrong as prose.
fck_apa_method <- function(m) {
  if (is.null(m) || !length(m)) return("--")
  switch(tolower(as.character(m)[1]),
         bonferroni = "Bonferroni", holm = "Holm", hochberg = "Hochberg",
         hommel = "Hommel", bh = , fdr = "Benjamini-Hochberg",
         by = "Benjamini-Yekutieli", tukey = "Tukey HSD",
         `games-howell` = "Games-Howell", none = "no",
         as.character(m)[1])
}

fck_md_table <- function(df, align = NULL) {
  if (is.null(df) || !nrow(df)) return(character(0))
  hdr <- names(df)
  if (is.null(align)) align <- c("l", rep("r", length(hdr) - 1))
  sep <- vapply(align, function(a) if (identical(a, "l")) ":---" else "---:",
                character(1))
  rows <- apply(df, 1, function(r) paste0("| ", paste(r, collapse = " | "), " |"))
  c(paste0("| ", paste(hdr, collapse = " | "), " |"),
    paste0("| ", paste(sep, collapse = " | "), " |"),
    unname(rows))
}


# ---- the report --------------------------------------------------------------
#
# Returns a character vector of Markdown lines. `values` and `input` are the
# app's own, and nothing is recomputed here, so the report cannot disagree with
# what was on screen.
fck_apa_report <- function(values, input, title = NULL) {
  L <- character(0)
  add   <- function(...) L <<- c(L, paste0(...))
  blank <- function() L <<- c(L, "")
  h     <- function(level, text) { blank(); add(strrep("#", level), " ", text); blank() }
  `%|%` <- function(a, b) if (is.null(a) || (length(a) == 1 && is.na(a))) b else a

  pkg <- function(p) tryCatch(as.character(utils::packageVersion(p)),
                              error = function(e) NA_character_)
  n_sub  <- if (!is.null(values$data)) nrow(values$data) else NA_integer_
  n_time <- if (!is.null(values$data)) ncol(values$data) else NA_integer_

  # ---------------------------------------------------------------- header --
  add("# ", title %|% "Functional data analysis: methods and results")
  blank()
  add("*Generated by F\\*CK on ", format(Sys.time(), "%Y-%m-%d %H:%M"), ".*")
  blank()
  add("> **How to use this document.** The *Statistical analysis* section is ",
      "written to be adapted into a Methods section and *Results* into Results, ",
      "in APA 7th-edition style. Every number comes from the analyses actually ",
      "run in this session; nothing is recomputed here, so the report cannot ",
      "disagree with what the app displayed. **Read the \"What these numbers do ",
      "not establish\" note that closes each Results subsection before ",
      "reporting anything** &mdash; several of the quantities this app produces ",
      "are descriptive and are commonly, and wrongly, reported as inferential.")
  blank()

  # ------------------------------------------------------------------ data --
  h(2, "Data")
  if (is.na(n_sub)) {
    add("No data had been imported when this report was generated.")
    return(L)
  }
  add("The analysis matrix contained ", n_sub, " curve",
      if (n_sub == 1) "" else "s", " observed at ", n_time, " time points",
      if (!is.null(values$time_labels))
        paste0(", from ", values$time_labels[1], " to ",
               values$time_labels[length(values$time_labels)]) else "", ".")

  if (!is.null(values$subject_ids)) {
    n_id <- length(unique(values$subject_ids))
    if (n_id < n_sub)
      add("These rows come from ", n_id, " distinct participants; ",
          n_sub - n_id, " row", if (n_sub - n_id == 1) " is" else "s are",
          " a repeated observation of a participant already present. A ",
          "between-groups test would treat those as independent, which is ",
          "anticonservative; use a repeated-measures design where that matters.")
    else
      add("Each row is a distinct participant (*n* = ", n_id, ").")
  }

  # AUDIT (P13.3). This was the report's ONLY statement about missing data, and
  # it is gated on values$fill_status, which the smoothing step creates. Cosinor
  # is a regression on the observations and needs no smoothing, so a perfectly
  # ordinary workflow -- import, fit cosinor, write the report -- produced a
  # document that never mentioned missingness at all. APA 7 asks for the amount
  # and the handling of missing data whether or not an imputation step ran, and
  # a reader cannot judge a curve analysis without knowing how much of it was
  # measured. The raw matrix is always available, so the count is always
  # reportable; the observed/interpolated/extrapolated split still requires the
  # smoother, because only it knows which side of a curve's own range a filled
  # cell fell on.
  fs <- values$fill_status
  if (is.null(fs) && !is.null(values$data)) {
    dm <- as.matrix(values$data)
    n_cell <- length(dm); n_miss <- sum(is.na(dm))
    if (n_miss == 0) {
      add("No values were missing: every one of the ", n_cell,
          " participant-by-time cells was measured.")
    } else {
      rows_aff <- sum(apply(dm, 1, function(r) any(is.na(r))))
      add("Of ", n_cell, " participant-by-time cells, ", n_miss, " (",
          fck_apa_num(100 * n_miss / n_cell, 1), "%) were missing, affecting ",
          rows_aff, " of ", nrow(dm), " curve", if (rows_aff == 1) "" else "s",
          ". No smoothing step has been run in this session, so these were ",
          "handled by whichever analysis read them: the cosinor fits drop ",
          "missing rows per participant and are unbiased under MCAR, while any ",
          "curve method requires a smoothed representation first.")
    }
  }
  if (!is.null(fs)) {
    n_cell <- length(fs)
    cnt <- function(k) sum(fs == k, na.rm = TRUE)
    add("Of ", n_cell, " cells, ", cnt("observed"), " (",
        fck_apa_num(100 * cnt("observed") / n_cell, 1), "%) were measured, ",
        cnt("interpolated"), " (", fck_apa_num(100 * cnt("interpolated") / n_cell, 1),
        "%) interpolated between observations and ", cnt("extrapolated"), " (",
        fck_apa_num(100 * cnt("extrapolated") / n_cell, 1),
        "%) extrapolated beyond a curve's first or last observation. ",
        "Extrapolated values are constrained by the smoother rather than ",
        "supported by data.")
  }

  if (!is.null(values$group_variables) && length(values$group_variables)) {
    gv <- values$group_variables
    rows <- do.call(rbind, lapply(names(gv), function(nm) {
      f <- droplevels(as.factor(gv[[nm]]))
      data.frame(Variable = nm, Levels = nlevels(f),
                 `n per level` = paste(sprintf("%s = %d", levels(f),
                                               as.integer(table(f))), collapse = ", "),
                 check.names = FALSE, stringsAsFactors = FALSE)
    }))
    blank(); add("**Grouping variables available.**"); blank()
    L <- c(L, fck_md_table(rows, c("l", "r", "l")))
  }

  # ------------------------------------------------ statistical analysis -----
  h(2, "Statistical analysis")
  add("Analyses were carried out in R ", as.character(getRversion()),
      " with the *fda* package (version ", pkg("fda") %|% "--", ")",
      if (!is.na(pkg("mgcv"))) paste0(", *mgcv* (", pkg("mgcv"), ")") else "",
      if (!is.na(pkg("minpack.lm"))) paste0(" and *minpack.lm* (", pkg("minpack.lm"), ")") else "",
      ". All package versions are listed under Reproducibility below and are ",
      "recorded in the exported analysis script.")
  blank()

  sm <- values$smooth_fit_metrics
  if (!is.null(sm)) {
    meth <- sm$method %|% "auto"
    add("**Smoothing.** Each series was represented in a ",
        if (isTRUE(input$is_cyclic)) "Fourier" else "B-spline", " basis of ",
        sm$n_basis %|% "--", " functions over ", sm$time_axis %|% "the column index",
        ". ",
        if (identical(meth, "none"))
          "No roughness penalty was applied: the basis interpolates the observations, so each curve passes through every measured value."
        else if (identical(meth, "auto"))
          # P12.1: this used to say "over the whole sample", which was not true
          # -- the search capped its objective at 60 subjects. The cap is gone,
          # and rather than assert a scope the sentence now STATES the count, so
          # a reader can check it against the reported n. They can still differ
          # legitimately: a subject with too few observed points cannot be
          # scored and is excluded from the objective.
          paste0("The smoothing parameter was chosen by minimising the mean ",
                 "generalised cross-validation score across ",
                 if (!is.null(sm$lambda_n_used) && is.finite(sm$lambda_n_used))
                   paste0("all ", sm$lambda_n_used, " ",
                          if (!is.null(sm$n_subjects) &&
                              is.finite(sm$n_subjects) &&
                              sm$n_subjects != sm$lambda_n_used)
                            paste0("of the ", sm$n_subjects, " ") else "",
                          "subjects with enough observed time points to be scored")
                 else "the sample",
                 ", giving lambda = ",
                 formatC(sm$lambda %|% NA, format = "e", digits = 2),
                 ". A single parameter was used for the whole sample rather ",
                 "than one per subject, so the curves remain comparable. ",
                 "(In *fda*, lambda = 0 is the unpenalised fit rather than an ",
                 "automatic selection; the value quoted here is the one the GCV ",
                 "search returned.)")
        else
          paste0("The smoothing parameter was set by hand to lambda = ",
                 formatC(sm$lambda %|% NA, format = "e", digits = 2), "."),
        " Missing time points were filled by the smoother rather than by ",
        "listwise deletion, so every curve contributes over the whole interval.")
    blank()
  }

  wr <- values$warping_results
  if (!is.null(wr)) {
    mlab <- switch(wr$method %|% "",
      linear_shift = "shift registration, a translation of the time axis",
      parametric   = sprintf("parametric registration in the %s family",
                             wr$family %|% "selected"),
      landmark     = paste0("landmark registration, a monotone piecewise-linear map through matched landmarks",
                            if (!is.null(wr$landmarks_used) && length(wr$landmarks_used))
                              sprintf(" placed by hand at %s",
                                      paste(fck_apa_num(wr$landmarks_used, 3), collapse = ", "))
                            else " detected automatically, by the same deterministic rule for every curve (the strongest turning point of that curve nearest each reference position taken from the sample mean)"),
      wr$method %|% "an unnamed method")
    add("**Curve registration.** Curves were registered by ", mlab, ". ")
    if (identical(wr$method, "linear_shift"))
      add("A shift is a translation, not an endpoint-preserving ",
          "reparameterisation. ",
          if (identical(wr$boundary, "periodic wrap"))
            "The design was treated as periodic, so the warp wraps and every curve stays within the observed domain."
          else
            "The design was not treated as periodic, so values beyond the observed domain were filled by constant extrapolation; the proportion extrapolated is reported with the diagnostics.")
    else if (identical(wr$method, "parametric"))
      add("The parameter search was restricted to the interval on which the ",
          "family is a strictly increasing bijection of the observed interval, ",
          "and widened where necessary to contain that family's identity, so ",
          "that \"this curve needs no registration\" was always an available ",
          "answer.")
    if (identical(wr$method, "landmark")) {
      # P13.4: state the rule, not only its exceptions. "0 curves were rejected"
      # is a result; silence is not, and a reader cannot tell the two apart.
      if (!is.null(wr$n_rejected) && wr$n_rejected > 0)
        add(" Landmarks that cross or duplicate cannot define a monotone time ",
            "warp; ", wr$n_rejected, " curve", if (wr$n_rejected == 1) "" else "s",
            " yielded such landmarks and ", if (wr$n_rejected == 1) "was" else "were",
            " left UNREGISTERED rather than registered with a fold, and ",
            if (wr$n_rejected == 1) "is" else "are",
            " therefore carried through the analyses in their original time.")
      else
        add(" Landmarks that cross or duplicate cannot define a monotone time ",
            "warp, and such curves are left unregistered rather than registered ",
            "with a fold; no curve required that here.")
    }
    # P11.4: the estimate has a resolution, and a paper is the last place a
    # difference below it should be read as a finding. The report states the
    # floor next to the method, in the Methods section, where a reader checks
    # whether a reported phase difference is resolvable at all.
    ng <- length(wr$time_points %|% numeric(0))
    if (identical(wr$method, "linear_shift") && ng > 1)
      add(" The shift was estimated on the ", ng, "-point analysis grid, so ",
          "the resolution of the estimate is one grid step, ",
          fck_apa_num(1 / (ng - 1), 4, bounded = TRUE),
          " of the observed interval; differences smaller than that are not ",
          "resolvable and are not interpreted.")
    blank()
  }

  fa <- values$fanova_results
  if (!is.null(fa)) {
    add("**Functional ANOVA.** ")
    # AUDIT (P13.4). The report described the test and never said WHICH CURVES it
    # was computed on. With a registration in the session that is the difference
    # between a claim about amplitude alone and a claim about amplitude and phase
    # together, and a reader has no way to tell the two apart from the output.
    # `used_warped` is what the analysis actually ran on, not what was requested.
    if (!is.null(values$warping_results) || isTRUE(fa$used_warped)) {
      if (isTRUE(fa$used_warped))
        add("This test was computed on the REGISTERED curves",
            if (!is.null(fa$warp_method_used))
              paste0(" (", fa$warp_method_used, " registration)") else "",
            ". Registration removes phase variation into the warping functions, ",
            "so a group difference found here is a difference in AMPLITUDE at ",
            "aligned time; it is not evidence that the groups differ in timing. ",
            "Differences in timing are the warping parameters, and are tested ",
            "separately if at all. ")
      else
        add("This test was computed on the UNREGISTERED curves, even though a ",
            "registration exists in this session, so a group difference here ",
            "may reflect a difference in timing, in amplitude, or in both. ")
    }
    if (identical(fa$design, "within"))
      add("A pointwise repeated-measures ANOVA was computed at each evaluation ",
          "point on the complete cases available there. The null distribution ",
          "was obtained by permuting condition labels *within each ",
          "participant*, one relabelling per participant per replicate applied ",
          "at every time point, so that whole trajectories rather than ",
          "individual values were exchanged (*B* = ", format(fa$n_permutations),
          "). Effect size is partial eta-squared, formed from the same sums of ",
          "squares as the *F* statistic.")
    else
      add("A pointwise one-way ANOVA was computed at each evaluation point, ",
          "with the null distribution obtained by permuting group labels ",
          "across whole curves (*B* = ", format(fa$n_permutations),
          "). Effect size is eta-squared.")
    add(" Pointwise *p* values were adjusted across time points by the ",
        "Benjamini-Hochberg false discovery rate and evaluated at alpha = ",
        fck_apa_num(fa$alpha %|% .05, 2, bounded = TRUE),
        ". Monte Carlo *p* values were computed as (1 + #{*T*\\* >= *T*}) / ",
        "(1 + *B*), so they can never be exactly zero; the smallest attainable ",
        "value here is ", fck_apa_pval(1 / ((fa$n_permutations %|% 0) + 1), 4), ".")
    blank()
  }

  pw <- values$pairwise_results
  if (!is.null(pw)) {
    m <- length(pw$pair_names)
    add("**Post-hoc comparisons.** All ", m, " pairwise comparisons were tested ",
        "by the same permutation procedure (*B* = ", format(pw$n_permutations),
        if (!is.null(pw$omnibus_permutations) &&
            !identical(as.integer(pw$omnibus_permutations), as.integer(pw$n_permutations)))
          paste0("; the omnibus test used ", format(pw$omnibus_permutations)) else "",
        "), with the global statistic given by the integrated squared ",
        "difference between the two mean curves. Familywise error was ",
        "controlled by the ", fck_apa_method(pw$correction_method), " correction.")
    afl <- if (identical(pw$correction_method, "none")) 1 / (pw$n_permutations + 1)
           else m / (pw$n_permutations + 1)
    add(" The smallest attainable *adjusted* *p* at this permutation count is ",
        fck_apa_pval(afl, 4),
        if (afl > (pw$alpha %|% .05))
          ", which exceeds the alpha level used: at this resolution no comparison could reach significance whatever the data showed."
        else ".")
    blank()
  }

  pa <- values$pca_anova
  if (!is.null(pa)) {
    add("**Component scores.** Group differences on the functional principal ",
        "component scores were tested component by component. Omnibus *p* ",
        "values were adjusted across components by the ",
        fck_apa_method(pa$across_pc_correction), " method, and pairwise comparisons within a ",
        "component were computed only where the across-component adjusted ",
        "omnibus *p* passed a gate of ",
        fck_apa_num(pa$posthoc_gate %|% .05, 2, bounded = TRUE),
        ", using the ", fck_apa_method(pa$posthoc_method), " correction.")
    blank()
  }

  hm <- values$harmonic_model
  if (!is.null(hm)) {
    add("**Cosinor (harmonic) regression.** Each participant's series was ",
        "fitted with a cosinor model of period ", format(hm$period), " and ",
        hm$n_harmonics, " harmonic", if (hm$n_harmonics == 1) "" else "s",
        if (!identical(hm$trend_type, "none"))
          paste0(", with a ", switch(hm$trend_type,
                 exp_sat = "saturating-exponential", linear = "linear",
                 log = "logarithmic", hm$trend_type), " trend term") else "",
        ", fitted to the ",
        if (isTRUE(hm$using_smoothed)) "smoothed curves" else "raw observations", ". ",
        if (isTRUE(hm$using_smoothed))
          "Fitting on smoothed rather than raw data removes independent noise and induces residual autocorrelation: *R*-squared is inflated, leave-one-out cross-validation is optimistic, and the zero-amplitude *F* test is anticonservative. Consider re-running on the raw series and reporting the difference."
        else
          "Cosinor is a regression on the observations and handles missing time points natively, so no smoothing was required.")
    # AUDIT (P12.2). Everything below was missing, and each item is something a
    # chronobiology reviewer expects to see. Cornelissen (2014) and the
    # CosinorOnline reporting convention ask for the period, the reference time
    # the acrophase is measured from, MESOR/amplitude/acrophase with intervals,
    # the zero-amplitude test, and the percentage rhythm. The app computes all
    # of them; the report simply did not print them.
    if (!is.null(hm$dv_name))
      add(" The dependent variable was ", hm$dv_name,
          if (!is.null(hm$dv_units)) paste0(" (", hm$dv_units, ")") else "", ".")
    add(" The period was FIXED at ", format(hm$period),
        " rather than estimated, so the acrophase and amplitude are conditional ",
        "on that choice.")
    if (!is.null(fck_cosinor_origin_label(hm$time_origin)))
      add(" Acrophase is reported in clock time, measured from ",
          fck_cosinor_origin_label(hm$time_origin), "; a phase has no meaning ",
          "without the origin it is measured from.")
    add(" Rhythm detection used the zero-amplitude *F* test, which compares the ",
        "full model against a refit containing the trend but no harmonics; it ",
        "tests whether the rhythmic amplitude is distinguishable from zero, ",
        "not whether the cosinor is the correct functional form.")
    add(" Joint confidence regions for amplitude and acrophase were obtained ",
        "from the error ellipse of the (cosine, sine) coefficient pair rather ",
        "than by treating the two as independent, following Bingham et al. ",
        "(1982); where the ellipse contains the origin the acrophase is not ",
        "identified at that level and is reported as such rather than as a ",
        "point estimate.")
    if (!is.null(hm$group_var_name)) {
      # AUDIT (P13.1). This was `pw <- ...`, which REBOUND the name the fANOVA
      # post-hoc object is held in (assigned once near the top of the function
      # and read much further down, in Results). The cosinor Methods paragraph
      # runs in between, so the post-hoc Results block was reading a cosinor
      # data frame: it printed "0 of 0 comparisons remained significant after
      # the -- correction at alpha = --" on a run whose Methods paragraph, four
      # lines above, correctly described 6 comparisons at B = 777. Introduced by
      # me at P12 and shipped. The name is scoped to the cosinor block now, and
      # tests/testthat/test-p13-corrections.R generates a report with BOTH kinds
      # of pairwise result present, which is the only arrangement that shows it.
      hp_pw <- values$hp_pairwise_results
      add(" Participants were grouped by ", hm$group_var_name, ".")
      if (!is.null(hp_pw) && nrow(hp_pw) > 0) {
        prm  <- fck_cosinor_param_label(values$hp_pairwise_param %|% "")
        circ <- grepl("acro", values$hp_pairwise_param %|% "")
        add(" Group differences in ", prm, " were tested pairwise on the ",
            "per-participant estimates",
            if (circ)
              paste0(", using the Watson-Williams *F* test for circular means, ",
                     "which assumes a von Mises distribution with a common and ",
                     "sufficiently high concentration in each group")
            else ", using Welch's *t* test, which does not assume equal variances",
            ". ",
            local({
              cr <- values$hp_pairwise_correction %|% "none"
              if (identical(cr, "none"))
                "No correction was applied for the number of comparisons, so the pairwise *p*-values are nominal."
              else paste0("*p*-values were adjusted for the ", nrow(hp_pw),
                          " comparisons by the ", fck_apa_method(cr), " method.")
            }))
      }
    }
    blank()
  }

  rg <- values$reg_model
  if (!is.null(rg)) {
    add("**Function-on-scalar regression.** ")
    if (identical(rg$inference, "analytic-t-fdr")) {
      add("A separate ordinary least-squares regression was fitted at each time ",
          "point. Standard errors and *p* values are analytical, referred to ",
          "*t* on ", format(rg$df_resid), " degrees of freedom and adjusted ",
          "across time points by the Benjamini-Hochberg false discovery rate.")
      if (!is.null(rg$n_boot))
        add(" A residual bootstrap (*B* = ", format(rg$n_boot), ") supplied a ",
            "percentile confidence interval for the coefficient curves. It is ",
            "not the test: the residual bootstrap resamples from the same ",
            "homoscedastic model the analytical standard error assumes, so it ",
            "estimates the same quantity with additional Monte Carlo variance. ",
            "The interval and the test may therefore disagree at the margin.")
    } else {
      add("A generalised additive model was fitted over the (participant, time) ",
          "long form, with a penalised spline in time and a by-smooth per ",
          "predictor, estimated by REML. Coefficient curves are prediction ",
          "contrasts; no pointwise standard error is propagated through them, ",
          "so no coefficient-curve inference is reported. Term-level tests are ",
          "available from the model summary.")
    }
    blank()
  }

  cl <- values$clustering_results
  if (!is.null(cl)) {
    add("**Functional clustering.** Curves were clustered by ",
        cl$method_label %|% cl$method, " into ", cl$k, " clusters",
        if (!identical(cl$standardize %|% "none", "none"))
          paste0(", after ", cl$standardize, " standardisation") else "",
        ". Partition quality is described by the average silhouette width and ",
        "the Calinski-Harabasz index. Neither is a test, and the number of ",
        "clusters was set by the analyst rather than estimated.")
    blank()
  }
  # ------------------------------------------------------------- results -----
  h(2, "Results")
  any_result <- FALSE
  tlab <- function(frac) {
    if (is.null(values$time_labels) || !is.finite(frac)) return(NULL)
    k <- length(values$time_labels)
    values$time_labels[max(1, min(k, round(frac * (k - 1) + 1)))]
  }

  if (!is.null(sm) && is.finite(sm$mean_r_squared %|% NA)) {
    any_result <- TRUE
    h(3, "Smoothing fit")
    add("The smoother reproduced the observed series closely, ",
        "*R*<sup>2</sup> = ", fck_apa_num(sm$mean_r_squared, 3, bounded = TRUE),
        if (is.finite(sm$sd_r_squared %|% NA))
          paste0(" (*SD* = ", fck_apa_num(sm$sd_r_squared, 3, bounded = TRUE), ")") else "",
        ", *RMSE* = ", fck_apa_num(sm$mean_rmse, 3),
        if (is.finite(sm$sd_rmse %|% NA))
          paste0(" (*SD* = ", fck_apa_num(sm$sd_rmse, 3), ")") else "", ".",
        if (is.finite(sm$mean_df %|% NA))
          paste0(" The fitted curves used ", fck_apa_num(sm$mean_df, 1),
                 " effective degrees of freedom on average.") else "")
    blank()
    add("*What these numbers do not establish.* A high smoothing ",
        "*R*<sup>2</sup> means the basis followed the observations, not that ",
        "the smoothing parameter was well chosen: an interpolating fit attains ",
        "*R*<sup>2</sup> = 1.00 by construction.")
  }

  if (!is.null(fa)) {
    any_result <- TRUE
    h(3, if (identical(fa$design, "within"))
           "Repeated-measures functional ANOVA" else "Functional ANOVA")
    ns <- sum(fa$sig_regions, na.rm = TRUE); nt <- length(fa$p_values_adjusted)
    add("Comparing ", fa$n_groups, " ",
        if (identical(fa$design, "within")) "conditions" else "groups",
        if (!is.null(fa$group_var)) paste0(" of ", fa$group_var) else "", " (",
        paste(sprintf("%s: *n* = %d", fa$groups, fa$group_sizes), collapse = "; "),
        "), the curves differed at ", ns, " of ", nt,
        " evaluation points after false-discovery-rate adjustment (",
        fck_apa_num(100 * ns / nt, 1), "% of the interval).")
    if ((fa$n_undefined %|% 0) > 0)
      add(" *F* was undefined at ", fa$n_undefined, " further point",
          if (fa$n_undefined == 1) "" else "s",
          " where there was no residual variation; those points are reported ",
          "as missing rather than as null results.")
    idx <- which.max(replace(fa$F_stat, !is.finite(fa$F_stat), -Inf))
    if (length(idx) == 1 && is.finite(fa$F_stat[idx])) {
      when <- tlab(fa$time_points[idx])
      add(" The largest difference occurred",
          if (!is.null(when)) paste0(" near ", when) else "", ", ",
          fck_apa_F(fa$F_stat[idx], fa$df_between, fa$df_within,
                    fa$p_values_adjusted[idx]), " (FDR-adjusted).")
    }
    add(" ", if (identical(fa$eta_squared_type, "partial"))
              "Partial eta-squared" else "Eta-squared",
        " averaged ", fck_apa_num(mean(fa$eta_squared, na.rm = TRUE), 3, bounded = TRUE),
        " across the interval (range ", fck_apa_range(fa$eta_squared, 3, TRUE), ").")
    if (is.finite(fa$L2_stat %|% NA))
      add(" The global *L*<sup>2</sup> statistic was ",
          fck_apa_num(fa$L2_stat, 3),
          if (is.finite(fa$p_value_L2 %|% NA))
            paste0(", ", fck_apa_p(fa$p_value_L2), ".")
          else ". This procedure computes no global *p* value for a repeated-measures design; the optional global test is that of Kurylo and Smaga (2023).")
    blank()
    add("*What these numbers do not establish.* The pointwise procedure ",
        "answers *where* the curves differ, not *whether* they differ overall. ",
        "The false-discovery-rate adjustment controls the expected proportion ",
        "of false positives among the flagged points, not the familywise error. ",
        if (identical(fa$eta_squared_type, "partial"))
          "Partial eta-squared removes between-participant variance from its denominator and is therefore not comparable with the classical eta-squared reported for between-groups designs."
        else "")
  }

  if (!is.null(pw)) {
    any_result <- TRUE
    h(3, "Post-hoc pairwise comparisons")
    rows <- do.call(rbind, lapply(pw$pair_names, function(nm) {
      r <- pw$results[[nm]]
      data.frame(Comparison = nm, `n₁` = r$n1, `n₂` = r$n2,
                 `L²` = fck_apa_num(r$L2_stat, 3),
                 `p` = fck_apa_pval(r$p_value_L2),
                 `p adj` = fck_apa_pval(r$p_value_L2_adjusted),
                 check.names = FALSE, stringsAsFactors = FALSE)
    }))
    L <- c(L, fck_md_table(rows)); blank()
    nsig <- sum(vapply(pw$results, function(x) isTRUE(x$sig_global), logical(1)))
    add(nsig, " of ", length(pw$pair_names), " comparisons remained significant ",
        "after the ", fck_apa_method(pw$correction_method), " correction at alpha = ",
        fck_apa_num(pw$alpha, 2, bounded = TRUE), ".")
    fl <- 1 / (pw$n_permutations + 1)
    praw <- vapply(pw$results, function(x) as.numeric(x$p_value_L2), numeric(1))
    n_at <- sum(is.finite(praw) & praw <= fl + 1e-12)
    if (n_at > 0)
      add(" ", n_at, " comparison", if (n_at == 1) "" else "s",
          " reached the smallest *p* value this permutation count can produce (",
          fck_apa_pval(fl, 4), "); ", if (n_at == 1) "its true value is" else "their true values are",
          " below that bound and this run cannot separate ",
          if (n_at == 1) "it" else "them",
          ". Identical *p* values in the table above are the resolution limit, ",
          "not a tie.")
    blank()
    add("*What these numbers do not establish.* A permutation *p* value is a ",
        "Monte Carlo estimate: re-running with a different seed moves it by ",
        "roughly its own standard error, and no permutation test can report a ",
        "*p* below 1/(*B* + 1).")
  }

  if (!is.null(values$pca_results)) {
    any_result <- TRUE
    pr <- values$pca_results
    nk <- fck_n_harmonics(pr)
    if (nk >= 1) {
      h(3, "Functional principal component analysis")
      vp <- pr$varprop[seq_len(nk)]
      add("The first ", nk, " functional principal component",
          if (nk == 1) "" else "s", " accounted for ",
          fck_apa_num(100 * sum(vp), 1), "% of the variance",
          if (!is.null(wr))
            " of the *registered* curves, so these components describe amplitude variation; the phase variation is carried by the warping parameters"
          else "", ".")
      rows <- data.frame(Component = paste0("PC", seq_len(nk)),
                         `Variance (%)` = fck_apa_num(100 * vp, 1),
                         `Cumulative (%)` = fck_apa_num(100 * cumsum(vp), 1),
                         check.names = FALSE, stringsAsFactors = FALSE)
      blank(); L <- c(L, fck_md_table(rows)); blank()
      add("*What these numbers do not establish.* The number of components ",
          "retained was chosen by the analyst. Variance explained is not ",
          "evidence that a component is interpretable or that it would ",
          "reappear in a new sample.")
    }
  }

  if (!is.null(pa) && !is.null(pa$k)) {
    any_result <- TRUE
    h(3, "Group differences on component scores")
    rows <- do.call(rbind, lapply(seq_len(pa$k), function(j) {
      o <- pa$omnibus[[j]]
      data.frame(Component = paste0("PC", j),
                 `F` = if (is.null(o)) "--" else fck_apa_num(o$F, 2)[1],
                 # fck_pc_omnibus() names these df1/df2, not df_between/df_within;
                 # sprintf() on a NULL yields character(0) and data.frame() then
                 # fails with "arguments imply differing number of rows".
                 `df` = if (is.null(o)) "--" else sprintf("%s, %s", o$df1, o$df2),
                 `p` = fck_apa_pval(pa$p_omnibus[j]),
                 `p adj` = fck_apa_pval(pa$p_omnibus_adj[j]),
                 `ω²` = if (is.null(o)) "--" else fck_apa_num(o$omega2, 3, bounded = TRUE)[1],
                 check.names = FALSE, stringsAsFactors = FALSE)
    }))
    L <- c(L, fck_md_table(rows)); blank()
    add("*p* values are adjusted across components by the ",
        fck_apa_method(pa$across_pc_correction), " method.")
    blank()
    add("*What these numbers do not establish.* A component that separates the ",
        "groups is not thereby the component that matters scientifically: the ",
        "components were extracted without reference to the grouping, and their ",
        "order reflects variance, not relevance.")
  }

  if (!is.null(hm) && !is.null(hm$individual_fits)) {
    any_result <- TRUE
    h(3, "Cosinor regression")
    fits <- hm$individual_fits
    ok <- vapply(fits, function(f) isTRUE(f$success), logical(1))
    g1 <- function(nm) vapply(fits[ok], function(f) {
      v <- f[[nm]]; if (length(v)) as.numeric(v)[1] else NA_real_ }, numeric(1))
    add(sum(ok), " of ", length(fits), " series were fitted successfully.")
    if (sum(ok) > 0) {
      add(" MESOR ", fck_apa_M(g1("mesor")), "; amplitude ",
          fck_apa_M(g1("amplitudes")), "; acrophase ",
          fck_apa_M(g1("acrophases_time")), " h.")
      pv <- g1("p_value"); nsig <- sum(is.finite(pv) & pv < .05)
      add(" A detectable rhythm (zero-amplitude *F* test, alpha = .05) was ",
          "present in ", nsig, " of ", sum(ok), " series (",
          fck_apa_num(100 * nsig / sum(ok), 1), "%). Model fit was ",
          "*R*<sup>2</sup> ", fck_apa_M(g1("r_squared"), 3, bounded = TRUE), ".")
    }

    # ---- the joint confidence regions the app already computes --------------
    # AUDIT (P12.2): fck_bingham_ci() builds the error ellipse for every
    # subject and the report never printed it. An amplitude without an interval
    # and an acrophase without one are not publication-grade, and the
    # "identified" count is the honest way to say how often the ellipse
    # contained the origin -- i.e. how often the phase is not estimable at all.
    bs <- hm$bingham_summary
    if (!is.null(bs) && length(bs) >= 1 && !is.null(bs[[1]])) {
      b1 <- bs[[1]]
      add(" Joint 95% confidence regions (first harmonic) had a median ",
          "half-width of ", fck_apa_num(b1$median_amp_halfwidth, 3),
          " for the amplitude",
          if (is.finite(b1$median_acro_halfwidth_h %|% NA_real_))
            paste0(" and ", fck_apa_num(b1$median_acro_halfwidth_h, 2),
                   " h for the acrophase") else "",
          ". The acrophase was identified at the 95% level in ",
          b1$n_identified, " of ", b1$n, " series",
          if (b1$n_identified < b1$n)
            "; in the remainder the confidence ellipse contained the origin, so no phase is estimable at that level"
          else "", ".")
    }
    blank()

    # ---- per-group description ----------------------------------------------
    gf <- hm$group_fits
    if (!is.null(hm$group_var_name) && !is.null(gf) && length(gf) >= 1) {
      gf <- gf[!vapply(gf, function(g) isTRUE(g$is_unassigned), logical(1))]
    }
    if (!is.null(hm$group_var_name) && length(gf) >= 1) {
      add("**By ", hm$group_var_name, ".**")
      blank()
      grows <- do.call(rbind, lapply(names(gf), function(nm) {
        g <- gf[[nm]]
        data.frame(
          Group      = nm,
          `n`        = g$n %|% NA,
          `MESOR`    = fck_apa_num(g$rhythm_adjusted_mean %|% g$intercept %|% NA_real_, 2),
          `Amplitude`= paste0(fck_apa_num(g$mean_amplitude %|% NA_real_, 2), " (",
                              fck_apa_num(g$sd_amplitude %|% NA_real_, 2), ")"),
          `Acrophase (h)` = fck_apa_num(g$mean_acrophase_time %|% NA_real_, 2),
          check.names = FALSE, stringsAsFactors = FALSE)
      }))
      L <- c(L, fck_md_table(grows)); blank()
      add("MESOR is the rhythm-adjusted mean. Amplitude is given as *M* (*SD*) ",
          "of the per-participant estimates. Acrophase is a group mean of a ",
          "circular quantity and is reported for description only; the test ",
          "below is the inferential statement.")
      blank()
    }

    # ---- pairwise group comparisons ------------------------------------------
    # P12.2: every comparison the user ran, not just the last one. A cosinor
    # paper reports MESOR, amplitude AND acrophase; the app kept only the most
    # recent, so the report showed one of the three.
    pw_all <- values$hp_pairwise_all
    if (is.null(pw_all) && !is.null(values$hp_pairwise_results))
      pw_all <- list(list(results = values$hp_pairwise_results,
                          param = values$hp_pairwise_param %|% "",
                          correction = values$hp_pairwise_correction %|% "none"))
    for (entry in pw_all) {
      # P13.1: named for the loop it belongs to, not `hp_pw` -- see the note in the
      # Methods block above for what reusing that name cost.
      hp_pw <- entry$results
      if (is.null(hp_pw) || !nrow(hp_pw)) next
      prm  <- fck_cosinor_param_label(entry$param %|% "")
      circ <- grepl("acro", entry$param %|% "")
      cr   <- entry$correction %|% "none"
      add("**Group differences in ", prm, ".**")
      blank()
      prows <- do.call(rbind, lapply(seq_len(nrow(hp_pw)), function(i) {
        r <- hp_pw[i, ]
        data.frame(
          Comparison = r$comparison,
          `Group 1`  = sprintf("%s (%.2f, n = %d)", r$group1, r$mean1, as.integer(r$n1)),
          `Group 2`  = sprintf("%s (%.2f, n = %d)", r$group2, r$mean2, as.integer(r$n2)),
          Statistic  = if (circ)
                         sprintf("*F*(1, %d) = %s", as.integer(r$df),
                                 fck_apa_num(r$t_stat, 2))
                       else
                         sprintf("*t*(%s) = %s", fck_apa_num(r$df, 1),
                                 fck_apa_num(r$t_stat, 2)),
          `p`        = fck_apa_pval(r$p_adjusted),
          # AUDIT (P12.2): for the circular comparison this column used to print
          # r$cohens_d, which in that branch holds the difference in MEAN
          # RESULTANT LENGTHS -- a difference in concentration. Watson-Williams
          # tests a difference in MEANS, so the reported "effect" answered a
          # different question from the test beside it, and on concentrated data
          # (the usual case for acrophase) it printed -0.00 for every row while
          # the means differed by hours. The effect for a difference of circular
          # means is the angular difference itself, which the app already
          # computes on the shortest arc and which is what a chronobiology paper
          # reports. Delta-r is still described in the note below, as what it is.
          `Effect`   = if (circ) paste0(fck_apa_num(r$mean_diff, 2), " h")
                       else fck_apa_num(r$cohens_d, 2),
          check.names = FALSE, stringsAsFactors = FALSE)
      }))
      L <- c(L, fck_md_table(prows)); blank()
      add("*p* is ", if (identical(cr, "none")) "unadjusted"
          else paste0("adjusted by the ", fck_apa_method(cr), " method"),
          " across the ", nrow(hp_pw), " comparisons of this parameter",
          if (length(pw_all) > 1)
            "; the correction was applied within each parameter, not across parameters, so the family being controlled is the set of pairwise contrasts for one parameter"
          else "",
          ". The effect column is ",
          if (circ)
            paste0("the difference between the group circular means, on the ",
                   "shortest arc, in hours -- the quantity the test above is ",
                   "about. (The app also reports the difference in mean ",
                   "resultant length on screen; that measures a difference in ",
                   "CONCENTRATION, not in mean phase, and is not an effect size ",
                   "for this test.)")
          else "Cohen's *d*.")
      blank()
    }

    add("*What these numbers do not establish.* Acrophase is a **circular** ",
        "quantity. The mean quoted above is arithmetic and is interpretable ",
        "only when the acrophases do not straddle the period boundary; use the ",
        "circular summary on the polar-plot tab for a directional mean and its ",
        "concentration.",
        if (isTRUE(hm$using_smoothed))
          " Because the fits were computed on smoothed curves, the *R*-squared and the proportion of detectable rhythms are both optimistic."
        else "")
    if (length(pw_all) > 0) {
      params_run <- vapply(pw_all, function(e) e$param %|% "", character(1))
      circ <- any(grepl("acro", params_run))
      # AUDIT (P12.2). Two limitations that a chronobiology reviewer WILL raise,
      # and which the app was silent about.
      add(" These comparisons are a two-stage procedure: a cosinor is fitted ",
          "per participant and the resulting point estimates are then compared ",
          "between groups. That is not the population-mean cosinor of Bingham ",
          "et al. (1982), and it treats each participant's estimate as if it ",
          "were measured without error, so a participant whose own rhythm is ",
          "poorly determined counts as much as one whose rhythm is precise.")
      if (circ)
        add(" The Watson-Williams test additionally assumes von Mises ",
            "distributions with a common, sufficiently high concentration; it ",
            "is unreliable when the phases are widely dispersed, and the ",
            "concentration is reported with the circular summary.")
      # Bingham's own caveat, checked rather than merely recited.
      if (any(grepl("amp", params_run))) {
        ac <- values$hp_acrophase_differs
        add(" Bingham et al. (1982) note that a difference in amplitude cannot ",
            "be interpreted when the groups also differ in acrophase, because ",
            "the amplitude is then estimated about different phases",
            if (isTRUE(ac))
              ": in these data the acrophase comparison WAS significant, so the amplitude result above should not be read as a difference in rhythm strength."
            else if (identical(ac, FALSE))
              ": the acrophase comparison on these data was not significant, so that caution does not apply here."
            else ". Run the same comparison on acrophase to check whether it applies here.")
      }
    }
  }

  if (!is.null(rg) && !is.null(rg$beta.hat)) {
    any_result <- TRUE
    h(3, "Function-on-scalar regression")
    bh <- rg$beta.hat
    if (identical(rg$inference, "analytic-t-fdr")) {
      rows <- do.call(rbind, lapply(seq_len(nrow(bh)), function(i) {
        pp <- rg$beta.p[i, ]
        data.frame(Coefficient = rownames(bh)[i],
                   `β range` = fck_apa_range(bh[i, ], 3),
                   `Points p adj < .05` = sprintf("%d of %d",
                                                  sum(pp < .05, na.rm = TRUE), length(pp)),
                   `Min p adj` = fck_apa_pval(min(pp, na.rm = TRUE)),
                   check.names = FALSE, stringsAsFactors = FALSE)
      }))
      L <- c(L, fck_md_table(rows)); blank()
      add("Each coefficient curve is summarised by the range of the estimate ",
          "over the interval and the number of time points at which the ",
          "FDR-adjusted *p* fell below .05.")
    } else {
      add("Coefficient curves were obtained as prediction contrasts from the ",
          "fitted generalised additive model; no pointwise inference on those ",
          "curves is available.")
    }
    if (!is.null(rg$r2_t))
      add(" Pointwise *R*<sup>2</sup> averaged ",
          fck_apa_num(mean(rg$r2_t, na.rm = TRUE), 3, bounded = TRUE),
          " (range ", fck_apa_range(rg$r2_t, 3, TRUE), ").")
    blank()
    # AUDIT (P13.2). This caveat was printed unconditionally, and it describes
    # the POINTWISE OLS fit: "each time point was fitted separately", "the
    # pointwise intervals are not simultaneous bands". Under the GAM branch both
    # halves are false -- the whole point of that branch is a penalised spline in
    # time, so the curves are not fitted independently, and the sentence above it
    # has just said that NO pointwise inference is available, so there are no
    # intervals for it to be talking about. A caveat that describes a different
    # estimator is not a small stylistic problem: it tells a reader their
    # coefficient curves are unsmoothed when they are penalised, and warns them
    # about the coverage of bands that do not exist.
    if (identical(rg$inference, "analytic-t-fdr")) {
      add("*What these numbers do not establish.* Each time point was fitted ",
          "separately, so the coefficient curves carry no smoothness constraint ",
          "and the pointwise intervals are not simultaneous bands: the ",
          "probability that the whole curve lies inside them is lower than the ",
          "nominal level at any single point. The false-discovery-rate ",
          "adjustment controls the expected proportion of false positives among ",
          "the flagged time points, not the familywise error.")
    } else {
      add("*What these numbers do not establish.* The coefficient curves are ",
          "prediction contrasts read off the fitted model, not estimated ",
          "parameters with a sampling distribution attached: no standard error ",
          "is propagated through them, so nothing in this subsection is a test ",
          "of where a predictor matters, and the curves should not be read as ",
          "if they carried confidence bands. The smoothness you see is imposed ",
          "by the penalty and its amount was chosen by REML from these same ",
          "data. For term-level inference use the model summary, whose ",
          "*p*-values are approximate because the smoothing parameters were ",
          "estimated rather than fixed.")
    }
  }

  if (!is.null(cl)) {
    any_result <- TRUE
    h(3, "Functional clustering")
    add("The ", cl$k, "-cluster solution had sizes ",
        paste(cl$cluster_sizes, collapse = ", "), ", accounting for ",
        fck_apa_num(100 * (cl$r_squared %|% NA), 1),
        "% of the total sum of squares.")
    if (is.finite(cl$silhouette_width %|% NA))
      add(" The average silhouette width was ",
          fck_apa_num(cl$silhouette_width, 3, bounded = TRUE),
          if (is.finite(cl$ch_index %|% NA))
            paste0(" and the Calinski-Harabasz index ",
                   fck_apa_num(cl$ch_index, 1)) else "", ".")
    blank()
    add("*What these numbers do not establish.* k-means partitions any data ",
        "set, including one with no cluster structure at all. Silhouette width ",
        "and the Calinski-Harabasz index describe the partition obtained; ",
        "neither tests whether clusters exist.")
  }

  if (!is.null(wr) && !is.null(wr$fit_statistics)) {
    any_result <- TRUE
    st <- wr$fit_statistics$summary
    h(3, "Registration diagnostics")
    g <- st$dispersion_reduction
    add("Registration reduced between-curve dispersion by ",
        if (!is.finite(g %|% NA)) "an unreported amount"
        else paste0(fck_apa_num(100 * g, 1), "%"),
        " (*V*<sub>pre</sub> = ", fck_apa_num(st$total_dispersion_pre, 3),
        ", *V*<sub>post</sub> = ", fck_apa_num(st$total_dispersion_post, 3),
        "). The time axis moved by ",
        fck_apa_num(st$mean_phase_displacement, 4),
        " of the interval on average.")
    if (is.finite(g %|% NA) && g > 0.15 &&
        is.finite(st$mean_phase_displacement %|% NA) &&
        st$mean_phase_displacement < 0.02)
      add(" **Caution.** A large reduction obtained from a near-identity warp ",
          "usually means the criterion absorbed *amplitude* differences rather ",
          "than aligning phase &mdash; near a peak, a small move in time changes ",
          "the value a great deal. Inspect the registered curves before ",
          "reporting this as evidence of phase variation.")
    blank()
    add("*What these numbers do not establish.* This is **not** an ",
        "amplitude/phase variance decomposition: *V*<sub>pre</sub> and ",
        "*V*<sub>post</sub> are two dispersions of the same curves, before and ",
        "after registration, and are not orthogonal components of a total. No ",
        "information criterion is reported, because no probability model for ",
        "the observed curves under a candidate registration is specified and ",
        "there is therefore no likelihood to penalise. Do not report the ",
        "reduction as \"variance explained by phase\".")
  }

  if (!any_result) {
    add("No analysis had been run when this report was generated. Run an ",
        "analysis and regenerate it.")
  }

  # -------------------------------------------------------- reproducibility --
  h(2, "Reproducibility")
  pkgs <- c("shiny", "fda", "mgcv", "minpack.lm", "cluster", "fda.usc",
            "rmfanova", "plotly", "ggplot2", "dplyr", "DT", "readxl")
  rows <- do.call(rbind, c(
    list(data.frame(Package = "R", Version = as.character(getRversion()),
                    stringsAsFactors = FALSE)),
    lapply(pkgs, function(p) {
      v <- pkg(p)
      if (is.na(v)) NULL else data.frame(Package = p, Version = v,
                                         stringsAsFactors = FALSE) })))
  L <- c(L, fck_md_table(rows, c("l", "r"))); blank()
  add("The **Download Analysis Code (R)** button on this tab writes a script ",
      "that reproduces every estimate above in a plain R session. It carries ",
      "the app's own estimators verbatim rather than a re-implementation, and ",
      "records these versions so a later run reports any that have changed.")
  blank()
  add("This project is **not** environment-pinned until `renv.lock` has been ",
      "generated on the analysis machine (`Rscript tools/renv_bootstrap.R`).")
  blank()
  add("---")
  blank()
  add("*Prepared with F\\*CK. Check every number against the app before ",
      "submission: this report is an aid to writing, not a substitute for ",
      "reading your own results.*")

  L
}


# ==============================================================================
# fck_apa_html(md_lines, title) -- a standalone HTML rendering
# ==============================================================================
# Deliberately small: headings, bold/italic, tables, blockquotes, rules,
# paragraphs. It exists so the report can be opened, read and printed to PDF
# without pandoc or rmarkdown, neither of which this app depends on. The APA
# table rules (top, below the header, bottom; no vertical lines) are in the
# stylesheet, so the printed output is close to journal style.
fck_apa_html <- function(md_lines, title = "Analysis report") {
  # Report text contains user-supplied strings: column names, factor levels, the
  # title. Escape BOTH angle brackets unconditionally, then restore only the
  # exact tag forms the report itself writes. An allow-list applied by NOT
  # escaping (the earlier approach) leaves ">" through and depends on the
  # negative lookahead being exhaustive; escaping first and restoring a fixed
  # literal set cannot be widened by anything in the input.
  keep <- c("sub", "sup", "strong", "em", "b", "i")
  esc <- function(x) {
    x <- gsub("&(?!(amp|lt|gt|mdash|nbsp);)", "&amp;", x, perl = TRUE)
    x <- gsub("<", "&lt;", x, fixed = TRUE)
    x <- gsub(">", "&gt;", x, fixed = TRUE)
    for (t in keep) {
      x <- gsub(paste0("&lt;", t, "&gt;"),  paste0("<", t, ">"),  x, fixed = TRUE)
      x <- gsub(paste0("&lt;/", t, "&gt;"), paste0("</", t, ">"), x, fixed = TRUE)
    }
    gsub("&lt;br&gt;", "<br>", x, fixed = TRUE)
  }
  # A literal asterisk is written "\\*" in the Markdown (the app is called
  # F*CK). Leaving it in place while the italic pass runs makes the emphasis
  # regex mis-pair its delimiters, so the title rendered as "F\\CK ... 19:26.*".
  # Park it out of the way first and restore it last.
  inline <- function(x) {
    x <- gsub("\\\\\\*", "\u0001", x)                     # literal * -> sentinel
    x <- gsub("\\*\\*(.+?)\\*\\*", "<strong>\\1</strong>", x)
    x <- gsub("\\*([^*]+?)\\*", "<em>\\1</em>", x)
    x <- gsub("`([^`]+)`", "<code>\\1</code>", x)
    gsub("\u0001", "*", x)                                # sentinel -> literal *
  }
  cells <- function(r) trimws(strsplit(gsub("^\\||\\|$", "", r), "\\|")[[1]])

  out <- character(0); i <- 1L; n <- length(md_lines)
  while (i <= n) {
    ln <- md_lines[i]
    if (!nzchar(trimws(ln))) { i <- i + 1L; next }
    if (grepl("^#{1,6} ", ln)) {
      lev <- nchar(sub("^(#+).*$", "\\1", ln))
      out <- c(out, sprintf("<h%d>%s</h%d>", lev,
                            inline(esc(sub("^#+ ", "", ln))), lev))
      i <- i + 1L; next
    }
    if (grepl("^-{3,}$", trimws(ln))) { out <- c(out, "<hr>"); i <- i + 1L; next }
    if (grepl("^> ", ln)) {
      blk <- character(0)
      while (i <= n && grepl("^> ", md_lines[i])) {
        blk <- c(blk, sub("^> ", "", md_lines[i])); i <- i + 1L }
      out <- c(out, "<blockquote><p>",
               inline(esc(paste(blk, collapse = " "))), "</p></blockquote>")
      next
    }
    if (grepl("^\\|", ln) && i < n && grepl("^\\|[ :|-]+\\|$", md_lines[i + 1L])) {
      hdr <- cells(ln); i <- i + 2L; body <- character(0)
      while (i <= n && grepl("^\\|", md_lines[i])) {
        body <- c(body, paste0("<tr><td>",
                  paste(inline(esc(cells(md_lines[i]))), collapse = "</td><td>"),
                  "</td></tr>")); i <- i + 1L }
      out <- c(out, "<table><thead><tr><th>",
               paste(inline(esc(hdr)), collapse = "</th><th>"),
               "</th></tr></thead><tbody>", body, "</tbody></table>")
      next
    }
    para <- character(0)
    while (i <= n && nzchar(trimws(md_lines[i])) &&
           !grepl("^(#{1,6} |\\||> |-{3,}$)", md_lines[i])) {
      para <- c(para, md_lines[i]); i <- i + 1L }
    out <- c(out, "<p>", inline(esc(paste(para, collapse = " "))), "</p>")
  }

  c('<!DOCTYPE html>', '<html lang="en"><head><meta charset="utf-8">',
    sprintf("<title>%s</title>", esc(title)),
    "<style>",
    "body{font-family:Georgia,'Times New Roman',serif;font-size:12pt;",
    "line-height:1.9;max-width:44em;margin:3em auto;padding:0 1.5em;color:#1a1a1a;}",
    "h1{font-size:1.6em;margin:0 0 .2em;line-height:1.3;}",
    "h2{font-size:1.22em;margin:2.4em 0 .4em;border-bottom:1px solid #d8d6d1;",
    "padding-bottom:.25em;}",
    "h3{font-size:1.04em;margin:1.8em 0 .3em;font-style:italic;}",
    "p{margin:.7em 0;text-align:justify;}",
    "table{border-collapse:collapse;margin:1.4em 0;font-size:.94em;width:100%;}",
    "th,td{padding:.42em .7em;border:0;}",
    "thead th{border-top:1px solid #1a1a1a;border-bottom:1px solid #1a1a1a;",
    "text-align:left;font-weight:600;}",
    "tbody tr:last-child td{border-bottom:1px solid #1a1a1a;}",
    "td:not(:first-child),th:not(:first-child){text-align:right;}",
    "blockquote{margin:1.6em 0;padding:.9em 1.2em;background:#f6f5f2;",
    "border-left:3px solid #7a7873;font-size:.95em;}",
    "blockquote p{text-align:left;}",
    "hr{border:0;border-top:1px solid #d8d6d1;margin:2.6em 0;}",
    "code{font-family:ui-monospace,Menlo,Consolas,monospace;font-size:.88em;",
    "background:#f2f1ee;padding:.1em .32em;border-radius:3px;}",
    "@media print{body{margin:0;max-width:none;font-size:11pt;}",
    "h2,h3{page-break-after:avoid;} table,blockquote{page-break-inside:avoid;}}",
    "</style></head><body>", out, "</body></html>")
}
