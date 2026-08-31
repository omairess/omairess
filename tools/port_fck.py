#!/usr/bin/env python3
"""
port_fck.py — one-time port script, kept for provenance.

Builds the machine-portable parts of the merged F*CK app out of the two
source apps that were merged:

    WaPaa1_3.R   "Functional Data Analysis Suite"      (fPCA / time-warped PCA,
                                                        functional ANOVA,
                                                        functional clustering)
    CIRCAREG.R   "Functional Regression Suite"         (FoSR, SoFR,
                                                        harmonic/cosinor
                                                        regression)

Everything that the two apps did NOT share — every analysis tab and every
analysis output — is carried across VERBATIM by line range, so the merged app
computes and prints exactly what the originals did.  Everything the two apps
DID share — data import, variable selection, smoothing, smoothing
diagnostics — was unified by hand; those files live in FCK/ and are NOT
written by this script.

The app in FCK/ is the source of truth from here on.  This script exists so
that any line of ported code can be traced back to the source app and line
range it came from; re-running it only rewrites the files listed in MANIFEST
and never touches the hand-written ones.
"""

import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
WAPAA = os.path.join(ROOT, "WaPaa1_3.R")
CIRCA = os.path.join(ROOT, "CIRCAREG.R")
OUT = os.path.join(ROOT, "FCK")


def read(path):
    with open(path, encoding="utf-8") as fh:
        return fh.read().split("\n")


def slice_lines(lines, first, last):
    """1-based inclusive line range."""
    return lines[first - 1:last]


# ---------------------------------------------------------------------------
# Renames applied to the ported CIRCAREG "Pairwise Comparisons" module.
#
# WaPaa also has a tab called "pairwise" (post-hoc tests after functional
# ANOVA) using input$run_pairwise / input$pairwise_correction.  CIRCAREG's
# pairwise tab compares COSINOR PARAMETERS between groups — a different
# analysis that happened to use the same input ids.  In the merged app the
# cosinor one is prefixed hp_ ("harmonic pairwise") so both tabs coexist.
# Nothing about either analysis changes; only the ids do.
# ---------------------------------------------------------------------------
HP_VALUE_RENAMES = [
    (r"values\$pairwise_results", "values$hp_pairwise_results"),
    (r"values\$pairwise_param", "values$hp_pairwise_param"),
    (r"values\$pairwise_correction", "values$hp_pairwise_correction"),
]

HP_ID_RENAMES = [
    ("export_pairwise_results", "hp_export_results"),
    ("export_pairwise_plot", "hp_export_plot"),
    ("pairwise_show_effect_size", "hp_show_effect_size"),
    ("pairwise_show_ci", "hp_show_ci"),
    ("pairwise_matrix_help", "hp_matrix_help"),
    ("pairwise_matrix", "hp_matrix"),
    ("pairwise_correction", "hp_correction"),
    ("pairwise_param", "hp_param"),
    ("pairwise_results", "hp_results"),
    ("pairwise_plot", "hp_plot"),
    ("run_pairwise", "hp_run"),
]


def rename_hp(text):
    # The tab itself also has to move: WaPaa's post-hoc tab already owns
    # tabName "pairwise", and two tabItems with the same tabName means only
    # one of them is ever reachable from the sidebar.
    text = text.replace('tabName = "pairwise"', 'tabName = "harm_pairwise"')
    for pat, rep in HP_VALUE_RENAMES:
        text = re.sub(pat, rep.replace("\\", "\\\\"), text)
    for old, new in HP_ID_RENAMES:
        text = re.sub(r"\b%s\b" % re.escape(old), new, text)
    return text


def patch(text, anchor, replacement, path, required=True):
    """Replace exactly one occurrence of `anchor`; fail loudly otherwise."""
    n = text.count(anchor)
    if n != 1:
        if not required and n == 0:
            return text
        raise SystemExit(
            "port_fck.py: expected exactly one occurrence of anchor in %s, "
            "found %d:\n%s" % (path, n, anchor))
    return text.replace(anchor, replacement)


# ---------------------------------------------------------------------------
# Surgical, documented edits to ported code.  Each one is anchored on an exact
# string from the source so that a change upstream turns into a hard error
# rather than a silent mis-port.
# ---------------------------------------------------------------------------

# The harmonic tab detects clock times from the column names itself.  The
# merged app ALSO detects them once, at import, into values$time_numeric
# (WaPaa's extractor).  Rather than have two detectors disagree, the harmonic
# time selector gains one EXTRA choice that reuses the shared vector.  The
# default is unchanged ("_index_"), so the original behaviour is untouched
# unless the user picks the new option.
HARMONIC_TIME_UI_ANCHOR = """      selectInput("harmonic_time_var", "Time Variable:", 
                  choices = c("Use column index (equally spaced)" = "_index_", 
                              "Specify times manually" = "_manual_",
                              numeric_vars),
                  selected = "_index_"),"""

HARMONIC_TIME_UI_NEW = """      selectInput("harmonic_time_var", "Time Variable:", 
                  choices = c("Use column index (equally spaced)" = "_index_", 
                              "Specify times manually" = "_manual_",
                              # MERGED APP: reuse the clock times detected once
                              # at import (values$time_numeric) instead of
                              # re-detecting them here.  Additive: the default
                              # is still "_index_".
                              "Use shared times detected at import" = "_shared_",
                              numeric_vars),
                  selected = "_index_"),
      conditionalPanel(
        condition = "input.harmonic_time_var == '_shared_'",
        if(!is.null(values$time_numeric) && length(values$time_numeric) == n_time) {
          div(style = "color: green; font-size: 0.9em;", icon("check-circle"),
              sprintf(" Using the %d time values detected at import: %s%s",
                      n_time, paste(head(values$time_numeric, 6), collapse = ", "),
                      if(n_time > 6) ", ..." else ""))
        } else {
          div(style = "color: #b00; font-size: 0.9em;", icon("exclamation-triangle"),
              " Import did not detect usable time values from the column names. Use 'Specify times manually'.")
        }
      ),"""

HARMONIC_TIME_SERVER_ANCHOR = """      if(input$harmonic_time_var == "_index_") {"""

HARMONIC_TIME_SERVER_NEW = """      if(input$harmonic_time_var == "_shared_") {
        # MERGED APP: the shared import step already extracted numeric clock
        # times from the column names into values$time_numeric.  Use them
        # verbatim so the harmonic tab and the fPCA/fANOVA/clustering tabs
        # all place the same column at the same time.
        if(is.null(values$time_numeric) || length(values$time_numeric) != n_time) {
          showNotification(
            "No shared time values available (import did not detect times from the column names). Use 'Specify times manually'.",
            type = "error", duration = 10)
          return()
        }
        time_vec <- as.numeric(values$time_numeric)
        original_times <- time_vec

      } else if(input$harmonic_time_var == "_index_") {"""

# CIRCAREG's export tab had a stub "Download Reproduction R Code" button whose
# handler wrote a single placeholder comment line.  WaPaa has a real 550-line
# code generator; the merged Export tab uses that one, and the stub is dropped
# rather than shipped as a working-looking button that produces nothing.
# CIRCAREG's FoSR coefficient export keeps working under a non-colliding id
# (WaPaa's export_scores_csv is the fPCA scores export).
# WaPaa's landmark plot branches on input$landmark_target and
# input$selected_subject, but its UI never creates either input.  Both are
# therefore NULL, and `NULL == "mean"` is logical(0), which makes `||` throw
# "invalid length zero argument" on R >= 4.3 — i.e. the landmark plot errors
# out.  This is a pre-existing bug in the source app, carried into the merge;
# reordering the test so the NULL check comes first fixes the crash and leaves
# behaviour identical if the inputs are ever supplied.
LANDMARK_GUARD_ANCHOR = """      if(input$landmark_target == "mean" || is.null(input$selected_subject)) {"""

LANDMARK_GUARD_NEW = """      # MERGED APP: NULL-safe ordering (see tools/port_fck.py). Neither input
      # is created by any UI, so both are NULL and the original test raised
      # "invalid length zero argument" instead of drawing the mean curve.
      if(is.null(input$selected_subject) || is.null(input$landmark_target) ||
         input$landmark_target == "mean") {"""

FOSR_EXPORT_ANCHOR = """  output$export_scores_csv <- downloadHandler("""
FOSR_EXPORT_NEW = """  output$export_fosr_coefs_csv <- downloadHandler("""


def header(target, sources):
    lines = [
        "# " + "=" * 74,
        "# %s" % target,
        "#",
        "# PORTED VERBATIM by tools/port_fck.py — do not hand-edit the ranges",
        "# below without updating that script's manifest.  Provenance:",
    ]
    for src, first, last, note in sources:
        lines.append("#   %s lines %d-%d%s" % (src, first, last,
                                               ("  (%s)" % note) if note else ""))
    lines += ["# " + "=" * 74, ""]
    return lines


# ---------------------------------------------------------------------------
# MANIFEST: every generated file, and exactly where its code comes from.
#   (relative path, [(source key, first line, last line, note), ...], transform)
# ---------------------------------------------------------------------------
MANIFEST = [
    # ---- UI: tabs carried across unchanged --------------------------------
    ("ui/30_diagnostics.R", "ui_tab_smooth_diag",
     [("W", 232, 348, "Smoothing Diagnostics tab", None)]),
    ("ui/40_settings.R", "ui_tab_settings",
     [("W", 351, 460, "fPCA / time-warped PCA settings", None)]),
    ("ui/41_results.R", "ui_tab_results",
     [("W", 463, 576, "Functional PCA results", None)]),
    ("ui/50_fanova.R", "ui_tab_fanova",
     [("W", 579, 694, "Functional ANOVA", None)]),
    ("ui/51_posthoc.R", "ui_tab_posthoc",
     [("W", 697, 789, "fANOVA post-hoc tests", None)]),
    ("ui/60_clustering.R", "ui_tab_clustering",
     [("W", 792, 1065, "Functional clustering", None)]),
    ("ui/70_fosr.R", "ui_tab_fosr",
     [("C", 275, 343, "Function-on-Scalar regression", None)]),
    ("ui/71_sofr.R", "ui_tab_sofr",
     [("C", 346, 421, "Scalar-on-Function regression", None)]),
    ("ui/72_harmonic.R", "ui_tab_harmonic",
     [("C", 424, 586, "Harmonic (cosinor) regression", None)]),
    ("ui/73_cosinor_pairwise.R", "ui_tab_cosinor_pairwise",
     [("C", 589, 654, "cosinor pairwise group tests; ids prefixed hp_", "hp")]),

    # ---- server: helpers ---------------------------------------------------
    ("server/01_helpers_time.R", None,
     [("W", 1152, 1382, "clock-time helpers used by every plot", None)]),
    ("server/02_helpers_gam.R", None,
     [("C", 702, 774, "GAM prediction helpers used by FoSR", None)]),

    # ---- server: shared views around the unified import / smoothing --------
    # The import and smoothing OBSERVERS are hand-merged (FCK/server/10_import.R
    # and 20_smoothing.R) because each is a union of the two apps' versions.
    # The read-only views around them had no real conflict — WaPaa's are strict
    # supersets — so they come across verbatim, plus CIRCAREG's compact
    # fit-metrics panel which WaPaa had no equivalent of.
    ("server/11_import_views.R", None,
     [("W", 1817, 1843, "recommended n_basis when the data change", None),
      ("W", 1902, 1967, "raw data plot with clock-time axis", None)]),
    ("server/21_smoothing_views.R", None,
     [("W", 1406, 1494, "smoothing fit statistics printout", None),
      ("C", 913, 930, "compact smoothing fit-metrics panel", None),
      ("W", 2862, 3016, "interactive smoothed-curve plot + curve selection", None)]),

    # ---- server: shared smoothing diagnostics ------------------------------
    # Both apps had this section; WaPaa's is a strict superset of CIRCAREG's
    # (it adds stratified CV folds and the GCV-vs-n-basis sweep) and its
    # reactiveValues names are the ones the merged app uses, so WaPaa's is
    # carried across whole and CIRCAREG's duplicate is dropped.
    ("server/30_diagnostics.R", None,
     [("W", 2231, 2860, "GAM REML, REML profile, CV, n-basis sweep", None)]),

    # ---- server: analysis families ----------------------------------------
    ("server/40_fpca.R", None,
     [("W", 3018, 4423, "group summary, fPCA/warping analysis + outputs", None),
      ("W", 6229, 6522, "warping / alignment / landmark plots", "landmark")]),
    ("server/50_fanova.R", None,
     [("W", 4425, 6227, "group UIs, fANOVA (between + repeated measures)", None),
      ("W", 6524, 6976, "post-hoc pairwise outputs", None)]),
    ("server/60_clustering.R", None,
     [("W", 6978, 8912, "functional clustering, optimisation, DCF, outputs", None)]),
    ("server/70_fosr.R", None,
     [("C", 1658, 2242, "Function-on-Scalar regression", None)]),
    ("server/71_sofr.R", None,
     [("C", 2244, 2877, "Scalar-on-Function regression", None)]),
    ("server/72_harmonic.R", None,
     [("C", 2879, 7190, "cosinor core, harmonic regression + outputs", "harmonic")]),
    ("server/73_cosinor_pairwise.R", None,
     [("C", 7264, 7859, "cosinor pairwise group tests; ids prefixed hp_", "hp")]),
    ("server/90_export.R", None,
     [("W", 8914, 9929, "all WaPaa exports + the reproducible-code generator", None),
      ("C", 7192, 7198, "FoSR coefficient export -> export_fosr_coefs_csv", "fosr_export"),
      ("C", 7209, 7262, "harmonic parameter + summary exports", None)]),
]


def build():
    src = {"W": read(WAPAA), "C": read(CIRCA)}
    src_name = {"W": "WaPaa1_3.R", "C": "CIRCAREG.R"}

    for rel, ui_var, pieces in MANIFEST:
        body = []
        for key, first, last, _note, transform in pieces:
            chunk = "\n".join(slice_lines(src[key], first, last))
            if transform == "hp":
                chunk = rename_hp(chunk)
            elif transform == "harmonic":
                chunk = patch(chunk, HARMONIC_TIME_UI_ANCHOR,
                              HARMONIC_TIME_UI_NEW, rel)
                chunk = patch(chunk, HARMONIC_TIME_SERVER_ANCHOR,
                              HARMONIC_TIME_SERVER_NEW, rel)
            elif transform == "landmark":
                chunk = patch(chunk, LANDMARK_GUARD_ANCHOR,
                              LANDMARK_GUARD_NEW, rel)
            elif transform == "fosr_export":
                chunk = patch(chunk, FOSR_EXPORT_ANCHOR, FOSR_EXPORT_NEW, rel)
            elif transform is not None:
                raise SystemExit("unknown transform %r" % transform)
            body.append(chunk)
        text = "\n\n".join(body)

        if ui_var:
            # a tabItem(...) slice ends in ")," inside the tabItems() list;
            # here it becomes a standalone assignment, so drop the comma.
            text = text.rstrip()
            if text.endswith(","):
                text = text[:-1]
            text = "%s <- %s" % (ui_var, text.lstrip())

        hdr = header(rel, [(src_name[k], a, b, n) for k, a, b, n, _t in pieces])
        path = os.path.join(OUT, rel)
        os.makedirs(os.path.dirname(path), exist_ok=True)
        with open(path, "w", encoding="utf-8") as fh:
            fh.write("\n".join(hdr) + text.rstrip() + "\n")
        print("wrote %-34s %6d lines" % (rel, text.count("\n") + 1))


if __name__ == "__main__":
    if not (os.path.exists(WAPAA) and os.path.exists(CIRCA)):
        sys.exit("port_fck.py: WaPaa1_3.R and CIRCAREG.R must be in %s" % ROOT)
    build()
