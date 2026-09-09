# ==============================================================================
# server/06_helpers_posthoc.R — what exactly are the post-hoc tests comparing?
#
# Pure logic, no Shiny: dance_posthoc_spec() takes plain lists so it can be tested
# without a session (tests/posthoc_source_test.R). The controls that drive it
# live in server/52_posthoc_source.R.
#
# THE BUG THIS FIXES
# ------------------
# The omnibus functional ANOVA runs on the variable chosen in
# `input$fanova_group_var` (via get_fanova_group_labels()), optionally
# restricted to a subset of its levels, on an fd object subset to match. The
# post-hoc tests then called perform_pairwise_comparisons() with
# `values$group_labels` — the PRIMARY grouping variable, i.e. whichever scalar
# variable happened to be selected first at import.
#
# So with two or more scalar variables selected, the omnibus tested one
# variable and the post-hoc tests silently tested another, on the full
# unfiltered curve set. Nothing in the output said so.
#
# WHAT REPLACES IT
# ----------------
# The post-hoc tests now take an explicit source:
#
#   "fanova" (default)  the labels, the curves and the design the omnibus
#                       actually used — read back from values$fanova_results,
#                       so they cannot drift apart by construction
#   "custom"            a variable and design the user picks here
#
# "custom" is allowed because a paired comparison on a variable the omnibus was
# never run on is a legitimate thing to want. It is not, however, a post-hoc
# test: it is a fresh family of comparisons, and the app says so rather than
# letting it borrow the omnibus's authority.
# ==============================================================================

# Resolve what to compare. Returns a list with $ok and, when ok, everything the
# two comparison functions need. Never guesses: if something does not line up it
# returns $ok = FALSE and a message that names the problem.
dance_posthoc_spec <- function(input, values) {
  bad <- function(msg) list(ok = FALSE, message = msg)

  res <- values$fanova_results
  if (is.null(res)) return(bad("Run the functional ANOVA first."))

  source_mode <- input$posthoc_source
  if (is.null(source_mode)) source_mode <- "fanova"

  n_curves <- function(fd) if (is.null(fd)) NA_integer_ else ncol(fd$coefs)

  # ---- follow the omnibus --------------------------------------------------
  if (identical(source_mode, "fanova")) {
    # The fd the omnibus used, so a level-subset in the omnibus carries over.
    fd <- res$fd_used
    if (is.null(fd)) fd <- values$fd_obj
    if (is.null(fd)) return(bad("No curves available — apply smoothing first."))

    if (identical(res$design, "within")) {
      if (is.null(res$subject_id) || is.null(res$rm_factor))
        return(bad("The functional ANOVA result has no subject ID / repeated-measures factor."))
      return(list(ok = TRUE, design = "within", fd = fd,
                  subject_id = res$subject_id, rm_factor = res$rm_factor,
                  matches_omnibus = TRUE,
                  description = sprintf(
                    "Within-subjects (paired), on '%s' — the same factor, curves and subjects as the omnibus test.",
                    res$group_var %||% "the repeated-measures factor")))
    }

    labs <- res$group_labels
    if (is.null(labs)) return(bad("The functional ANOVA result did not record its group labels; re-run it."))
    if (!is.na(n_curves(fd)) && length(labs) != n_curves(fd))
      return(bad(sprintf("Group labels (%d) and curves (%d) do not line up; re-run the functional ANOVA.",
                         length(labs), n_curves(fd))))
    return(list(ok = TRUE, design = "between", fd = fd,
                group_labels = droplevels(as.factor(labs)),
                matches_omnibus = TRUE,
                description = sprintf(
                  "Between-subjects, on '%s' (%s) — the same variable, levels and curves as the omnibus test.",
                  res$group_var %||% "the fANOVA variable",
                  paste(levels(droplevels(as.factor(labs))), collapse = ", "))))
  }

  # ---- a variable chosen here ----------------------------------------------
  fd <- values$fd_obj
  if (is.null(fd)) return(bad("No curves available — apply smoothing first."))
  design <- input$posthoc_design
  if (is.null(design)) design <- "between"

  if (identical(design, "within")) {
    sv <- input$posthoc_subject_var; rv <- input$posthoc_rm_var
    if (is.null(sv) || is.null(rv) || !nzchar(sv) || !nzchar(rv))
      return(bad("Choose both a subject-ID variable and a repeated-measures factor."))
    if (is.null(values$covariates) || !all(c(sv, rv) %in% names(values$covariates)))
      return(bad("Those variables are not among the scalar variables selected at import."))
    subj <- values$covariates[[sv]]; rmf <- droplevels(as.factor(values$covariates[[rv]]))
    if (!is.na(n_curves(fd)) && length(subj) != n_curves(fd))
      return(bad(sprintf("The scalar variables have %d rows but there are %d curves.",
                         length(subj), n_curves(fd))))
    if (nlevels(rmf) < 2)
      return(bad(sprintf("'%s' has only %d level — nothing to compare.", rv, nlevels(rmf))))
    # a paired test needs subjects measured at more than one level
    per_subject <- tapply(as.character(rmf), subj, function(x) length(unique(x)))
    if (all(per_subject < 2))
      return(bad(sprintf("No subject appears at more than one level of '%s', so nothing can be paired.", rv)))
    return(list(ok = TRUE, design = "within", fd = fd,
                subject_id = subj, rm_factor = rmf,
                matches_omnibus = FALSE,
                description = sprintf(
                  "Within-subjects (paired), on '%s' with subjects from '%s' (%d levels; %d of %d subjects appear at more than one level).",
                  rv, sv, nlevels(rmf), sum(per_subject >= 2), length(per_subject))))
  }

  gv <- input$posthoc_group_var
  if (is.null(gv) || !nzchar(gv)) return(bad("Choose a grouping variable."))
  if (is.null(values$group_variables) || !(gv %in% names(values$group_variables)))
    return(bad("That variable is not among the scalar variables selected at import."))
  labs <- droplevels(as.factor(values$group_variables[[gv]]))
  if (!is.na(n_curves(fd)) && length(labs) != n_curves(fd))
    return(bad(sprintf("'%s' has %d values but there are %d curves.", gv, length(labs), n_curves(fd))))
  if (nlevels(labs) < 2)
    return(bad(sprintf("'%s' has only %d level — nothing to compare.", gv, nlevels(labs))))
  list(ok = TRUE, design = "between", fd = fd, group_labels = labs,
       matches_omnibus = FALSE,
       description = sprintf("Between-subjects, on '%s' (%s).",
                             gv, paste(levels(labs), collapse = ", ")))
}

# ==============================================================================
# Repeated-measures columns, aligned to the curves
#
# WaPaa's repeated-measures path reads its subject-ID and factor columns from
# values$uploaded_data — the raw imported frame, BEFORE the import step drops
# all-NA rows (and, in this app, rows below the minimum-measured-points
# threshold). If any row was dropped, those vectors are longer than the curve
# set and every subject is paired with the wrong curve. Nothing errors: the
# lengths only have to be plausible for the test to run and return numbers.
#
# values$covariates is rebuilt and row-filtered alongside values$data, so it is
# aligned by construction. Prefer it; fall back to the raw frame only when the
# column was never selected as a scalar variable, and then check the length
# rather than trusting it.
# ==============================================================================
dance_rm_column <- function(values, varname, n_expected = NULL) {
  if (is.null(varname) || !nzchar(varname)) return(NULL)
  if (is.null(n_expected) && !is.null(values$data)) n_expected <- nrow(values$data)

  if (!is.null(values$covariates) && varname %in% names(values$covariates)) {
    v <- values$covariates[[varname]]
    if (is.null(n_expected) || length(v) == n_expected) return(v)
  }
  if (!is.null(values$uploaded_data) && varname %in% names(values$uploaded_data)) {
    v <- values$uploaded_data[[varname]]
    if (!is.null(n_expected) && length(v) != n_expected) {
      stop(sprintf(
        paste("'%s' has %d values but there are %d curves. Rows were dropped at import",
              "(all-missing rows, or the minimum-measured-points threshold), so this",
              "column no longer lines up. Select it as a scalar variable on the Data",
              "Import tab and it will be filtered with the rest."),
        varname, length(v), n_expected))
    }
    return(v)
  }
  stop(sprintf("Variable '%s' not found.", varname))
}

# ==============================================================================
# SHARED PERMUTATION ARITHMETIC (R1, R7)
# ==============================================================================
# AUDIT (P20). The two pairwise-comparison kernels in server/50_fanova.R used
# TWO DIFFERENT p-value conventions in the same result object:
#
#   global  (1 + #{T* >= T}) / (1 + B)      -- the add-one Monte Carlo estimator
#   pointwise      #{T* >= T}  / B          -- the plain proportion
#
# The plain proportion is not an estimator of a p-value under a randomisation
# null: it can be exactly zero, which asserts an impossible event, and with the
# app's default B it is badly anti-conservative in the tail. With B = 49 it can
# report p = 0 at a point where the observed statistic was simply the largest of
# 50 exchangeable values -- the correct answer there is 1/50 = .02. So NO valid
# raw pointwise p can be below 1/(B+1), and any value under it is an artefact.
#
# Both kernels now call this, so the two conventions cannot drift apart again.
# `exceedances` counts permuted statistics at least as extreme as the observed
# one; `B` is the number of permutations DRAWN, not the number that came back
# finite -- a non-finite permuted statistic is a failure to evaluate the
# statistic, and dropping it from the denominator silently shrinks the reference
# set and inflates significance.
# `B` is normally one number. It is allowed to be a vector of the same length as
# `exceedances` for the pointwise fANOVA kernels, where a time point at which the
# statistic was undefined for some permutations has a SMALLER reference set than
# the number drawn -- the P4.7/P5.7 rule that such a permutation contributes no
# draw at that point rather than a silent zero. Anything else is a mistake and
# is caught here rather than recycled into wrong p-values.
dance_perm_p <- function(exceedances, B) {
  stopifnot(length(B) == 1L || length(B) == length(exceedances),
            all(is.finite(B)), all(B >= 0),
            all(is.finite(exceedances)), all(exceedances >= 0),
            all(exceedances <= B))
  (1 + exceedances) / (1 + B)
}

# The smallest p either kernel can return, for labelling the resolution floor.
dance_perm_p_floor <- function(B) 1 / (1 + B)

# ------------------------------------------------------------------------------
# STUDENTISING WITHOUT DESTROYING THE EXTREME CASE (R7)
# ------------------------------------------------------------------------------
# Both kernels formed t = mean / se and then wrote
#
#     t_stat[!is.finite(t_stat)] <- 0
#
# which silently converts the MOST extreme possible result into the least. Twelve
# participants whose paired difference is identical and non-zero have zero
# dispersion, so se = 0 and t = +Inf: perfect, unanimous separation. Coerced to
# zero it is reported as t = 0, d = 0, and -- because every permuted statistic is
# then also compared against zero -- p = 1. The one dataset where the answer is
# unambiguous is the one dataset the test cannot see.
#
# The two zero-denominator cases are not the same thing and must not be treated
# the same way:
#
#   0 / 0   no difference AND no dispersion. There is nothing to detect and the
#           statistic is genuinely undefined; 0 is the right value, and it
#           agrees with every permutation, so p -> 1 honestly.
#
#   c / 0   a non-zero difference with no dispersion. The statistic is +/-Inf.
#           This is a real value, it compares correctly against finite permuted
#           statistics (abs(Inf) >= abs(finite) is TRUE), and under sign-flipping
#           only the all-same-sign flips reproduce it -- so the p-value comes out
#           at the resolution floor, which is correct.
#
# NA in either argument is missing data rather than a degenerate denominator, and
# stays NA so the caller can see it.
dance_studentise <- function(numerator, denominator) {
  out <- numerator / denominator
  # NaN arises only from 0/0 (and from Inf/Inf, which cannot occur here because
  # a variance is finite); NA propagates from missing input and is left alone.
  out[is.nan(out)] <- 0
  out
}
