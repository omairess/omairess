# ==============================================================================
# server/42_fpca_anova.R — group comparisons on fPCA / warped-fPCA scores
#
# NEW MODULE (2026-09-03). The fPCA tab produced an n x k score matrix and
# stopped there: there was no way to ask whether the groups differ on any
# component. This adds a per-component ANOVA with post-hocs.
#
# The arithmetic is in server/09_helpers_pcanova.R (pure base R, testable
# without a session); this file is the controls, the reactive and the output.
#
# Design decisions worth knowing about, all of them visible in the output:
#
#   TWO multiplicity families. Pairwise-within-a-component and omnibus-across-
#   components are separate families and both are corrected, with independent
#   controls. Correcting only the pairwise tests is the common error and the
#   less important of the two: it is the k omnibus tests that manufacture a
#   "significant component" by chance.
#
#   The post-hoc gate runs on the ACROSS-COMPONENT-ADJUSTED omnibus p, not the
#   raw one, so a component that only looks significant before that correction
#   does not get a table of pairwise tests lending it credibility.
#
#   Welch by default when the variances differ. Brown-Forsythe decides, per
#   component, unless the user overrides. With n = 654 against n = 59 the
#   equal-variance F is not a safe default.
#
#   Eigenvalue separation is reported next to every component, because a group
#   difference on a component whose eigenvalue is nearly tied with its
#   neighbour's is a difference on an arbitrary rotation, not on that mode.
#
#   Repeated curves are detected and flagged. Scores are one row per CURVE; if
#   a participant contributes several, the between-groups test treats
#   correlated rows as independent and is anticonservative.
# ==============================================================================

# ---- controls ---------------------------------------------------------------
output$pca_anova_controls <- renderUI({
  if (is.null(values$pca_results)) {
    return(helpText("Run the fPCA (or time-warped fPCA) first."))
  }
  nc <- if (!is.null(values$pca_results$scores)) ncol(values$pca_results$scores) else 0
  if (nc < 1) return(helpText("The PCA produced no scores."))

  gv <- values$selected_group_vars
  if (is.null(gv) || !length(gv)) {
    return(helpText(HTML("No grouping variable is available. Select one or more
                          scalar variables at import to compare groups on the
                          component scores.")))
  }

  tagList(
    selectInput("pca_anova_group_var", "Group variable:",
                choices = gv, selected = gv[1]),
    sliderInput("pca_anova_ncomp", "Components to test:",
                min = 1, max = nc, value = min(nc, 4), step = 1),
    helpText(HTML("Every component tested enters the across-component
                   multiplicity family, so testing all of them costs power on
                   each. Test the ones you would interpret.")),
    hr(),
    selectInput("pca_anova_across", "Correction ACROSS components (the omnibus family):",
                choices = FCK_PC_CORRECTIONS[!FCK_PC_CORRECTIONS %in% c("tukey", "games-howell")],
                selected = "holm"),
    helpText(HTML("<b>This is the correction that matters most.</b> Testing k
                   components is k chances to find a difference. The scores are
                   orthogonal in sample but the tests are not independent -- all
                   k depend on the same estimated eigenfunctions -- so
                   Bonferroni/Holm are conservative here and BH assumes a
                   positive dependence that holds only approximately.")),
    hr(),
    selectInput("pca_anova_posthoc", "Correction WITHIN a component (the pairwise family):",
                choices = FCK_PC_CORRECTIONS, selected = "holm"),
    helpText(HTML("Tukey HSD assumes equal variances and a common n; with
                   unbalanced groups of very different spread prefer
                   <b>Games-Howell</b>, which uses Welch degrees of freedom per
                   pair. Both control the family rate by construction, so for
                   those two the reported p is already adjusted.")),
    numericInput("pca_anova_gate",
                 "Only run post-hocs where the adjusted omnibus p is below:",
                 value = 0.05, min = 0.0001, max = 1, step = 0.01),
    hr(),
    radioButtons("pca_anova_omnibus", "Omnibus test:",
                 choices = c("Automatic (Welch when variances differ)" = "auto",
                             "Fisher F (assumes equal variances)" = "fisher",
                             "Welch F (does not)" = "welch"),
                 selected = "auto"),
    numericInput("pca_anova_conf", "Confidence level for the intervals:",
                 value = 0.95, min = 0.5, max = 0.999, step = 0.01),
    hr(),
    actionButton("run_pca_anova", "Run component ANOVA",
                 class = "btn-primary", width = "100%")
  )
})

# ---- the analysis -----------------------------------------------------------
observeEvent(input$run_pca_anova, {
  req(values$pca_results)
  pca <- values$pca_results
  if (is.null(pca$scores)) {
    showNotification("The PCA produced no scores to test.", type = "error")
    return()
  }

  gvar <- input$pca_anova_group_var
  if (is.null(gvar) || !nzchar(gvar) ||
      is.null(values$group_variables) || !(gvar %in% names(values$group_variables))) {
    showNotification("Select a grouping variable.", type = "error")
    return()
  }
  labels_all <- values$group_variables[[gvar]]

  S <- as.matrix(pca$scores)
  nc <- min(as.integer(input$pca_anova_ncomp %||% 4), ncol(S))
  S <- S[, seq_len(nc), drop = FALSE]

  # The scores have one row per CURVE that entered the PCA. The grouping
  # variable has one row per imported subject. If those two have drifted apart
  # -- which happens whenever curves were filtered between import and the PCA --
  # a positional join silently pairs the wrong label with the wrong curve, and
  # nothing downstream would reveal it. Refuse rather than guess.
  if (length(labels_all) != nrow(S)) {
    showNotification(
      sprintf("Cannot align: the PCA has %d curves but '%s' has %d values. The score rows and the group labels are not the same set of subjects, so any join here would be positional guesswork.",
              nrow(S), gvar, length(labels_all)),
      type = "error", duration = NULL)
    return()
  }
  g <- droplevels(as.factor(labels_all))
  if (nlevels(g) < 2) {
    showNotification("The grouping variable has fewer than two levels present.",
                     type = "error")
    return()
  }

  use_welch <- switch(input$pca_anova_omnibus %||% "auto",
                      "auto" = NA, "fisher" = FALSE, "welch" = TRUE, NA)

  withProgress(message = "Testing components...", value = 0, {
    res <- fck_pc_anova_all(
      scores = S, g = g,
      eigenvalues = pca$values,
      varprop = pca$varprop,
      posthoc_method = input$pca_anova_posthoc %||% "holm",
      across_pc_correction = input$pca_anova_across %||% "holm",
      posthoc_gate = suppressWarnings(as.numeric(input$pca_anova_gate %||% 0.05)),
      conf = suppressWarnings(as.numeric(input$pca_anova_conf %||% 0.95)),
      subject_ids = values$subject_ids,
      use_welch = use_welch)
    incProgress(1)
  })

  if (is.null(res)) {
    showNotification("No component could be tested.", type = "error")
    return()
  }
  res$group_var <- gvar
  res$warped <- !is.null(values$warping_results)
  res$warping_method <- if (res$warped) values$warping_results$method else NULL

  # AUDIT: registration splits a curve into phase and amplitude. The scores are
  # the amplitude half; the warping functions are the phase half, and stopping
  # at the scores discards exactly the half that answers "do the groups differ
  # in TIMING". The warping parameters get the same machinery, as their own
  # multiplicity family -- pooling them with the components would let one
  # question borrow the other's correction.
  res$warp <- NULL
  if (res$warped) {
    wp <- fck_warp_params(values$warping_results)
    if (!is.null(wp) && nrow(wp) == length(g)) {
      res$warp <- fck_pc_anova_all(
        scores = as.matrix(wp), g = g,
        posthoc_method = input$pca_anova_posthoc %||% "holm",
        across_pc_correction = input$pca_anova_across %||% "holm",
        posthoc_gate = suppressWarnings(as.numeric(input$pca_anova_gate %||% 0.05)),
        conf = suppressWarnings(as.numeric(input$pca_anova_conf %||% 0.95)),
        subject_ids = values$subject_ids, use_welch = use_welch)
      if (!is.null(res$warp)) res$warp$names <- colnames(wp)
    } else if (!is.null(wp)) {
      res$warp_note <- sprintf(
        "The warping parameters have %d rows but %d curves entered the PCA, so they could not be aligned. No phase comparison is shown rather than a guessed one.",
        nrow(wp), length(g))
    } else {
      res$warp_note <- "The warping produced no parameter that varies between curves (every curve was left where it was), so there is nothing to compare."
    }
  }

  values$pca_anova <- res

  n_sig <- sum(res$p_omnibus_adj < res$posthoc_gate, na.rm = TRUE)
  showNotification(sprintf("Component ANOVA complete: %d of %d component(s) differ across '%s' after the across-component correction.",
                           n_sig, res$k, gvar),
                   type = "message", duration = 8)
})

# ---- the report -------------------------------------------------------------
output$pca_anova_results <- renderPrint({
  res <- values$pca_anova
  if (is.null(res)) {
    cat("Run the fPCA, choose a grouping variable, then press 'Run component ANOVA'.\n\n")
    cat("What this does: a one-way test of the group means on EACH principal\n")
    cat("component's scores, with post-hoc pairwise comparisons.\n\n")
    cat("Two multiplicity families are corrected separately -- the pairwise tests\n")
    cat("within a component, and the omnibus tests across components. The second\n")
    cat("is the one that manufactures a false 'significant component', and it is\n")
    cat("the one most analyses forget.\n")
    return(invisible(NULL))
  }

  hdr <- function(x) cat("\n", strrep("=", 74), "\n", x, "\n", strrep("=", 74), "\n", sep = "")

  cat("=== Group comparison on ", if (res$warped) "time-warped " else "", "fPCA scores ===\n\n", sep = "")
  cat("Grouping variable: ", res$group_var, "\n", sep = "")
  cat("Groups: ", paste(res$levels, collapse = ", "), "\n", sep = "")
  if (res$warped)
    cat("Warping method: ", res$warping_method, "\n",
        "  Scores are from the REGISTERED curves, so amplitude variation is being\n",
        "  compared with phase variation removed. A group difference that vanishes\n",
        "  after warping was a timing difference, not a magnitude one.\n", sep = "")
  cat("Components tested: ", res$k, "\n", sep = "")

  # ---- independence -------------------------------------------------------
  ic <- res$independence
  if (isTRUE(ic$checked) && !isTRUE(ic$ok)) {
    cat("\n! REPEATED CURVES DETECTED: ", ic$n_rows, " score rows from ", ic$n_subjects,
        " participants\n", sep = "")
    cat("  (", ic$n_repeated, " rows are repeat observations; up to ", ic$max_per_subject,
        " per participant).\n", sep = "")
    cat("  These tests treat every row as an independent observation. They are\n")
    cat("  therefore ANTICONSERVATIVE: the p-values below are too small by an amount\n")
    cat("  that grows with the within-participant correlation. Either average the\n")
    cat("  scores within participant before testing, or fit a mixed model with a\n")
    cat("  random intercept per participant. This is not corrected automatically\n")
    cat("  because which of those is right depends on your design.\n")
  }

  # ---- the omnibus table --------------------------------------------------
  hdr("Omnibus tests, one per component")
  cat("Across-component correction: ", res$across_pc_correction,
      "  (family size = ", res$k, ")\n\n", sep = "")
  cat(sprintf("%-5s %8s %10s %12s %10s %10s %9s  %s\n",
              "PC", "var%", "test", "statistic", "p", "p adj", "omega2", "separation"))
  cat(strrep("-", 92), "\n")
  for (j in seq_len(res$k)) {
    o <- res$omnibus[[j]]
    if (is.null(o)) { cat(sprintf("%-5s  (not testable)\n", paste0("PC", j))); next }
    het <- !is.null(res$variance[[j]]) && isTRUE(res$variance[[j]]$p < 0.05)
    used_welch <- is.finite(o$welch_p) && isTRUE(all.equal(res$p_omnibus[j], o$welch_p))
    sep <- if (!is.null(res$separation) && j <= nrow(res$separation))
      res$separation$stability[j] else ""
    cat(sprintf("%-5s %7s%% %10s %12s %10s %10s %9s  %s\n",
                paste0("PC", j),
                if (!is.null(res$varprop) && j <= length(res$varprop)) fmt1(100 * res$varprop[j]) else "-",
                if (used_welch) "Welch F" else "Fisher F",
                if (used_welch) sprintf("F(%s,%s)=%s", fmt1(o$welch_df1), fmt1(o$welch_df2), fmt2(o$welch_F))
                else sprintf("F(%d,%d)=%s", o$df1, o$df2, fmt2(o$F)),
                format.pval(res$p_omnibus[j], digits = 3, eps = 1e-16),
                format.pval(res$p_omnibus_adj[j], digits = 3, eps = 1e-16),
                fmt3(o$omega2), sep))
  }
  cat("\nomega2 is the less optimistic effect size (eta2 is reported per component\n")
  cat("below). With large n a tiny omega2 can carry a very small p: read both.\n")

  # ---- assumptions --------------------------------------------------------
  hdr("Assumption checks")
  for (j in seq_len(res$k)) {
    bf <- res$variance[[j]]; nm <- res$normality[[j]]
    if (is.null(bf) && is.null(nm)) next
    cat(sprintf("PC%d\n", j))
    if (!is.null(bf))
      cat(sprintf("  Equal variances (Brown-Forsythe): F(%d,%d) = %s, p = %s  -> %s\n",
                  bf$df1, bf$df2, fmt2(bf$F),
                  format.pval(bf$p, digits = 3, eps = 1e-16),
                  if (bf$p < 0.05) "variances DIFFER; Welch used" else "no evidence against"))
    if (!is.null(nm)) {
      cat(sprintf("  Residual skew = %s, excess kurtosis = %s", fmt2(nm$skew), fmt2(nm$excess_kurtosis)))
      if (isTRUE(nm$shapiro_applicable) && is.finite(nm$shapiro_p))
        cat(sprintf(", Shapiro-Wilk W = %s, p = %s", fmt3(nm$shapiro_W),
                    format.pval(nm$shapiro_p, digits = 3, eps = 1e-16)))
      cat("\n")
      if (isTRUE(nm$severe))
        cat("    ! Severely non-normal (|skew| > 2 or |excess kurtosis| > 7). The F test\n",
            "      is not reliable here; use the Kruskal-Wallis line reported below.\n", sep = "")
      else if (isTRUE(nm$shapiro_applicable) && is.finite(nm$shapiro_p) && nm$shapiro_p < 0.05)
        cat("    Shapiro rejects, but skew and kurtosis are mild: at this n Shapiro-Wilk\n",
            "    detects departures far too small to trouble the F test. Proceed.\n", sep = "")
    }
  }

  # ---- eigenvalue separation ---------------------------------------------
  if (!is.null(res$separation)) {
    hdr("Eigenvalue separation")
    cat("A component whose eigenvalue is nearly tied with its neighbour's has a\n")
    cat("barely-determined eigenfunction: the split between the two is an arbitrary\n")
    cat("rotation. A group difference there belongs to the PAIR, not to the component.\n\n")
    sp <- res$separation
    cat(sprintf("%-5s %14s %12s   %s\n", "PC", "eigenvalue", "rel. gap", "verdict"))
    for (i in seq_len(min(nrow(sp), res$k)))
      cat(sprintf("%-5s %14s %12s   %s\n", paste0("PC", sp$component[i]),
                  fmt4(sp$eigenvalue[i]), fmt3(sp$rel_gap[i]), sp$stability[i]))
  }

  # ---- per-component detail ----------------------------------------------
  for (j in seq_len(res$k)) {
    o <- res$omnibus[[j]]
    if (is.null(o)) next
    hdr(sprintf("PC%d in detail", j))

    cat("Group means on the PC", j, " scores:\n", sep = "")
    for (lv in o$levels)
      cat(sprintf("  %-16s %10s  (SD %s, n = %d)\n", lv,
                  fmt3(o$means[lv]), fmt3(o$sds[lv]), o$ns[lv]))
    cat("\n  A PC and its negation are the same component, so the SIGN of these means\n")
    cat("  is only interpretable against the plotted loading for PC", j, ".\n", sep = "")

    cat(sprintf("\n  Fisher F(%d,%d) = %s, p = %s\n", o$df1, o$df2, fmt2(o$F),
                format.pval(o$p, digits = 3, eps = 1e-16)))
    if (is.finite(o$welch_F))
      cat(sprintf("  Welch  F(%s,%s) = %s, p = %s\n", fmt1(o$welch_df1), fmt1(o$welch_df2),
                  fmt2(o$welch_F), format.pval(o$welch_p, digits = 3, eps = 1e-16)))
    if (is.finite(o$kw_chisq))
      cat(sprintf("  Kruskal-Wallis chi2(%d) = %s, p = %s   [rank-based, no normality assumption]\n",
                  o$kw_df, fmt2(o$kw_chisq), format.pval(o$kw_p, digits = 3, eps = 1e-16)))
    cat(sprintf("  eta2 = %s, omega2 = %s, epsilon2 = %s\n",
                fmt3(o$eta2), fmt3(o$omega2), fmt3(o$epsilon2)))

    ph <- res$posthoc[[j]]
    if (is.null(ph)) {
      cat(sprintf("\n  No post-hoc tests: the across-component adjusted omnibus p (%s) is above\n",
                  format.pval(res$p_omnibus_adj[j], digits = 3, eps = 1e-16)))
      cat(sprintf("  the gate of %s. Running them anyway would be testing pairs inside a\n",
                  fmt3(res$posthoc_gate)))
      cat("  family the omnibus did not open.\n")
      next
    }

    cat(sprintf("\n  Post-hoc pairwise (%s%s):\n", ph$method[1],
                if (isTRUE(ph$already_adjusted[1])) ", family-wise by construction" else ""))
    cat(sprintf("  %-14s %-14s %10s %20s %9s %10s %10s %9s\n",
                "group A", "group B", "diff", "CI", "stat", "p", "p adj", "Hedges g"))
    cat("  ", strrep("-", 102), "\n", sep = "")
    for (i in seq_len(nrow(ph))) {
      r <- ph[i, ]
      cat(sprintf("  %-14s %-14s %10s %20s %9s %10s %10s %9s\n",
                  r$a, r$b, fmt3(r$diff),
                  sprintf("[%s, %s]", fmt3(r$ci_lo), fmt3(r$ci_hi)),
                  sprintf("%s=%s", r$stat_label, fmt2(r$stat)),
                  format.pval(r$p_raw, digits = 3, eps = 1e-16),
                  if (isTRUE(r$already_adjusted)) "(same)" else format.pval(r$p_adj, digits = 3, eps = 1e-16),
                  fmt2(r$hedges_g)))
    }
    cat("\n  The interval is on the score scale, at the confidence level you set. A\n")
    cat("  difference whose interval spans zero is not distinguishable from none,\n")
    cat("  whatever the adjusted p says.\n")
  }

  # ---- the phase half ------------------------------------------------------
  if (!is.null(res$warp_note)) {
    hdr("Warping parameters")
    cat(res$warp_note, "\n")
  } else if (!is.null(res$warp)) {
    hdr("Warping parameters (the PHASE half of the registration)")
    cat("The component scores above describe AMPLITUDE variation, because they\n")
    cat("were computed on the registered curves. These are what registration took\n")
    cat("OUT -- how far each curve had to be moved to align it. If the groups\n")
    cat("differ in timing rather than in magnitude, this is where it shows.\n")
    cat("Corrected as their own family (", res$warp$k,
        " test(s), ", res$warp$across_pc_correction,
        "), not pooled with the components: 'do groups differ in phase' is a\n",
        "different question from 'do they differ in the k-th mode of amplitude'.\n", sep = "")

    cat(sprintf("\n%-52s %11s %12s %9s\n", "parameter", "p", "p adj", "omega2"))
    cat(strrep("-", 88), "\n")
    for (j in seq_len(res$warp$k)) {
      o <- res$warp$omnibus[[j]]
      nmj <- res$warp$names[j]
      lbl <- FCK_WARP_LABELS[[nmj]] %||% nmj
      if (is.null(o)) { cat(sprintf("%-52s  (not testable)\n", lbl)); next }
      cat(sprintf("%-52s %11s %12s %9s\n", lbl,
                  format.pval(res$warp$p_omnibus[j], digits = 3, eps = 1e-16),
                  format.pval(res$warp$p_omnibus_adj[j], digits = 3, eps = 1e-16),
                  fmt3(o$omega2)))
      for (lv in o$levels)
        cat(sprintf("    %-16s mean %10s (SD %9s, n = %d)\n", lv,
                    fmt4(o$means[lv]), fmt4(o$sds[lv]), o$ns[lv]))
      ph <- res$warp$posthoc[[j]]
      if (!is.null(ph)) {
        cat(sprintf("    %-14s %-14s %11s %22s %11s %9s\n",
                    "group A", "group B", "diff", "CI", "p adj", "Hedges g"))
        for (i in seq_len(nrow(ph))) {
          r <- ph[i, ]
          cat(sprintf("    %-14s %-14s %11s %22s %11s %9s\n", r$a, r$b, fmt4(r$diff),
                      sprintf("[%s, %s]", fmt4(r$ci_lo), fmt4(r$ci_hi)),
                      format.pval(if (isTRUE(r$already_adjusted)) r$p_raw else r$p_adj,
                                  digits = 3, eps = 1e-16),
                      fmt2(r$hedges_g)))
        }
      }
    }
    cat("\n  A shift is in the units of the registration's time axis (the app warps\n")
    cat("  on a 0-1 range, so 0.01 is 1% of the recording). Warp amplitude is the\n")
    cat("  RMS distance of h(t) from the identity and is defined for any method,\n")
    cat("  including nonlinear warps where a single shift does not exist.\n")
  }

  # ---- what has been corrected for ---------------------------------------
  hdr("Multiplicity, in full")
  n_ph <- sum(vapply(res$posthoc, function(z) if (is.null(z)) 0L else nrow(z), integer(1)))
  cat(sprintf("Family 1 -- omnibus across components: %d test(s), corrected by %s.\n",
              res$k, res$across_pc_correction))
  cat(sprintf("Family 2 -- pairwise within a component: %d test(s) in total, corrected by %s\n",
              n_ph, res$posthoc_method))
  cat("            (separately within each component, which is the standard\n")
  cat("             convention; correcting Family 2 across components as well would\n")
  cat("             be defensible and more conservative still).\n")
  cat("\nThe two families are NOT jointly corrected. The post-hoc gate is what links\n")
  cat("them: pairs are only examined inside a component the corrected omnibus\n")
  cat("already opened, which is the usual protected-test logic.\n")
  cat("\nDependence: the k score columns are orthogonal in sample, but the k tests\n")
  cat("share the same estimated eigenfunctions and the same curves. Holm and\n")
  cat("Bonferroni remain valid (they need no independence) but are conservative;\n")
  cat("BH assumes positive regression dependence, which is plausible here but not\n")
  cat("guaranteed; BY is valid under any dependence and is the safe FDR choice.\n")

  invisible(NULL)
})

# ---- a plot of the score distributions --------------------------------------
output$pca_anova_plot <- renderPlotly({
  res <- values$pca_anova
  validate(need(!is.null(res), "Run the component ANOVA to see the score distributions."))
  req(values$pca_results, values$group_variables)

  S <- as.matrix(values$pca_results$scores)
  g <- droplevels(as.factor(values$group_variables[[res$group_var]]))
  validate(need(length(g) == nrow(S), "Score rows and group labels are not aligned."))

  k <- res$k
  df <- do.call(rbind, lapply(seq_len(k), function(j) data.frame(
    pc = factor(sprintf("PC%d (%s)", j,
                        if (!is.null(res$varprop)) paste0(fmt1(100 * res$varprop[j]), "%") else "-"),
                levels = sprintf("PC%d (%s)", seq_len(k),
                                 if (!is.null(res$varprop)) paste0(fmt1(100 * res$varprop[seq_len(k)]), "%") else "-")),
    score = S[, j], grp = g, stringsAsFactors = FALSE)))

  # Categorical slots of the reference palette, in fixed order, never cycled.
  cols <- c("#2a78d6", "#eb6834", "#1baf7a", "#eda100",
            "#e87ba4", "#008300", "#4a3aa7", "#e34948")
  lv <- levels(g)

  p <- plot_ly()
  for (i in seq_along(lv)) {
    d <- df[df$grp == lv[i], ]
    p <- add_trace(p, type = "box", y = d$score, x = d$pc,
                   name = lv[i], legendgroup = lv[i],
                   marker = list(color = cols[((i - 1) %% length(cols)) + 1], size = 3),
                   line = list(color = cols[((i - 1) %% length(cols)) + 1]),
                   fillcolor = paste0(cols[((i - 1) %% length(cols)) + 1], "33"),
                   boxpoints = "outliers")
  }
  # mark the components the corrected omnibus called different
  sig <- which(is.finite(res$p_omnibus_adj) & res$p_omnibus_adj < res$posthoc_gate)
  ann <- lapply(sig, function(j) list(
    x = levels(df$pc)[j], y = max(S[, j], na.rm = TRUE),
    text = sprintf("p_adj = %s", format.pval(res$p_omnibus_adj[j], digits = 2, eps = 1e-16)),
    showarrow = FALSE, yshift = 14,
    font = list(size = 11, color = "#52514e")))

  plotly::layout(p, boxmode = "group",
                 title = sprintf("Component scores by %s%s", res$group_var,
                                 if (res$warped) " (time-warped)" else ""),
                 xaxis = list(title = ""),
                 yaxis = list(title = "component score"),
                 annotations = ann,
                 legend = list(orientation = "h", y = -0.12))
})

# ---- CSV of everything ------------------------------------------------------
output$download_pca_anova <- downloadHandler(
  filename = function() sprintf("fpca_component_anova_%s.csv", Sys.Date()),
  content = function(file) {
    res <- values$pca_anova
    if (is.null(res)) { utils::write.csv(data.frame(), file, row.names = FALSE); return() }
    omni <- do.call(rbind, lapply(seq_len(res$k), function(j) {
      o <- res$omnibus[[j]]; if (is.null(o)) return(NULL)
      data.frame(component = j,
                 var_explained_pct = if (!is.null(res$varprop)) 100 * res$varprop[j] else NA,
                 F_fisher = o$F, df1 = o$df1, df2 = o$df2, p_fisher = o$p,
                 F_welch = o$welch_F, p_welch = o$welch_p,
                 KW_chisq = o$kw_chisq, p_KW = o$kw_p,
                 p_used = res$p_omnibus[j], p_across_adjusted = res$p_omnibus_adj[j],
                 eta2 = o$eta2, omega2 = o$omega2, epsilon2 = o$epsilon2,
                 stringsAsFactors = FALSE)
    }))
    ph <- do.call(rbind, lapply(seq_len(res$k), function(j) {
      z <- res$posthoc[[j]]; if (is.null(z)) return(NULL)
      cbind(component = j, z)
    }))
    utils::write.csv(
      list(omnibus = omni, posthoc = ph)$omnibus, file, row.names = FALSE)
    if (!is.null(ph)) {
      utils::write.table("", file, append = TRUE, col.names = FALSE, row.names = FALSE)
      utils::write.table(ph, file, append = TRUE, sep = ",",
                         row.names = FALSE, col.names = TRUE)
    }
  }
)
