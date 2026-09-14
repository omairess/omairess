# ==============================================================================
# server/02c_helpers_design.R — reading the Study Design inputs
# ==============================================================================
# One tiny pure file, because several modules need this and none of them should
# have to load the trajectory layer (and lme4 with it) to ask which factor the
# user is grouping by. Sourced early; called at runtime, so order never matters
# for the app -- it matters for the tests, which source an explicit list.
# ==============================================================================

if (!exists("%||%", mode = "function")) `%||%` <- function(a, b) if (is.null(a)) b else a

# ------------------------------------------------------------------------------
# THE ONE GROUPING VARIABLE, DERIVED FROM THE STUDY DESIGN
# ------------------------------------------------------------------------------
# P21 phase 4 replaced a single "Group Variable" control with the Study Design
# panel. Everything written against the old input -- the fitted-curve plot, the
# polar plot, the parameter histograms, the pairwise module, the density rings --
# needs one grouping factor, and this is where that factor is decided, once.
#
# It is a PURE function of `input` rather than a reactive because several of its
# callers are pure functions that take `input` as an argument and are tested in
# isolation. A first attempt wired those call sites to a reactive and broke six
# tests immediately: a reactive cannot be called outside a reactive context, and
# a blanket rename across a file that mixes observers with pure helpers cannot
# tell the two apart. The reactive wrapper in the server body calls this.
#
# A between-participant factor wins when there is one; otherwise the
# within-participant factor, since these displays split on whatever single
# factor the design offers. "_none_" is the legacy sentinel every caller already
# tests for, so returning it keeps their existing guards correct.
dance_group_var_from_design <- function(input) {
  pick <- function(x) if (!is.null(x) && nzchar(x) && !identical(x, "_none_")) x else NULL
  d <- input$harmonic_design %||% "between"
  out <- switch(d,
    between = pick(input$harmonic_between_var),
    within  = pick(input$harmonic_within_var),
    mixed   = pick(input$harmonic_between_var2) %||% pick(input$harmonic_within_var2),
    NULL)
  out %||% "_none_"
}
