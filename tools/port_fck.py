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
# merged app ALSO parses them once, at import, into values$time_clock (see
# FCK/server/03_helpers_clock.R).  Rather than have two detectors disagree, the
# harmonic time selector gains one EXTRA choice that reuses the shared vector.
# The default is unchanged ("_index_"), so the original behaviour is untouched
# unless the user picks the new option.
#
# NOTE: this deliberately uses values$time_clock, NOT values$time_numeric.
# WaPaa's extract_time_values() returns 1:n_time without reading the column
# names at all, so time_numeric is a column counter; feeding it to a cosinor
# fit would silently claim unequally-spaced measurements were hourly.
HARMONIC_TIME_UI_ANCHOR = """      selectInput("harmonic_time_var", "Time Variable:", 
                  choices = c("Use column index (equally spaced)" = "_index_", 
                              "Specify times manually" = "_manual_",
                              numeric_vars),
                  selected = "_index_"),"""

HARMONIC_TIME_UI_NEW = """      selectInput("harmonic_time_var", "Time Variable:", 
                  choices = c("Use column index (equally spaced)" = "_index_", 
                              "Specify times manually" = "_manual_",
                              # MERGED APP: reuse the clock times parsed once
                              # at import (values$time_clock) instead of
                              # re-detecting them here.  Additive: the default
                              # is still "_index_".
                              "Use shared clock times parsed at import" = "_shared_",
                              numeric_vars),
                  selected = "_index_"),
      conditionalPanel(
        condition = "input.harmonic_time_var == '_shared_'",
        if(!is.null(values$time_clock) && length(values$time_clock) == n_time) {
          div(style = "color: green; font-size: 0.9em;", icon("check-circle"),
              sprintf(" Using the %d clock times parsed at import: %s%s",
                      n_time, paste(head(values$time_clock, 6), collapse = ", "),
                      if(n_time > 6) ", ..." else ""))
        } else {
          div(style = "color: #b00; font-size: 0.9em;", icon("exclamation-triangle"),
              " Import could not parse clock times from the column names. Use 'Specify times manually'.")
        }
      ),"""

HARMONIC_TIME_SERVER_ANCHOR = """      if(input$harmonic_time_var == "_index_") {"""

HARMONIC_TIME_SERVER_NEW = """      if(input$harmonic_time_var == "_shared_") {
        # MERGED APP: the shared import step already parsed real clock hours
        # from the column names into values$time_clock (03_helpers_clock.R).
        # Use those, so this tab and the manual-entry route agree.
        if(is.null(values$time_clock) || length(values$time_clock) != n_time) {
          showNotification(
            "No shared clock times available: the column names did not yield hours in [0, 24). Use 'Specify times manually'.",
            type = "error", duration = 10)
          return()
        }
        time_vec <- as.numeric(values$time_clock)
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

# The smoothing fit summary reports the method, basis count and lambda but not
# which time axis was used.  Both source apps only ever had one (the column
# index), so there was nothing to report; the merged app lets the user smooth
# against real clock times instead, and a fit summary that does not say which
# axis produced it is a summary you cannot act on.
TIME_AXIS_ANCHOR = """        cat("Method:", method_label, "\\n")
      }"""

TIME_AXIS_NEW = """        cat("Method:", method_label, "\\n")
      }

      # MERGED APP: which time axis these numbers came from.
      if(!is.null(metrics$time_axis)) {
        cat("Time axis:", metrics$time_axis, "\\n")
      }"""

# Both the fPCA tab and the fANOVA tab build values$fd_obj themselves when the
# user goes straight to an analysis without visiting the smoothing tab.  Each
# picked min(20, n_time - 2) basis functions and said nothing about it, so on
# 24 hourly columns an analysis the user believed was running on raw data was
# running on a smooth.  Both now go through fck_ensure_fd_obj()
# (FCK/server/04_helpers_fd.R): one rule, an interpolating basis, and a notice
# saying which representation is in use.
# The cross-validation diagnostic fixes the basis at min(20, n_time - 2) while
# CIRCAREG's version used the user's n_basis.  The CV curve exists to recommend
# a smoothing factor for the smoothing the user is actually going to run, so a
# basis count unrelated to theirs makes the recommendation mis-targeted.  This
# follows n_basis (CIRCAREG's behaviour), still capped at n_time - 2.
# ---------------------------------------------------------------------------
# C1: the reproducible-code export.  WaPaa's generator covered its own family
# only; CIRCAREG had no generator at all (its button wrote a placeholder).  To
# emit code for the cosinor / FoSR / SoFR fits, the generator needs the
# settings each fit ACTUALLY ran with -- reading input$... at export time would
# emit whatever the widgets happen to say now, which is not what produced the
# stored model.  So each fit records its own settings alongside itself.
# ---------------------------------------------------------------------------
FOSR_SETTINGS_ANCHOR = """      values$reg_model <- fit
      updateSelectInput(session, "reg_color_var", choices = input$reg_predictors)"""

FOSR_SETTINGS_NEW = """      # MERGED APP: what this fit actually used, for the code export.
      fit$fck_settings <- list(
        predictors   = input$reg_predictors,
        method       = input$reg_method,
        use_bootstrap = isTRUE(input$use_bootstrap),
        n_boot       = if(isTRUE(input$use_bootstrap)) input$n_boot else NULL,
        using_smoothed = !is.null(values$smooth_data),
        n_subjects   = nrow(Y),
        n_time       = ncol(Y))

      values$reg_model <- fit
      updateSelectInput(session, "reg_color_var", choices = input$reg_predictors)"""

SOFR_SETTINGS_ANCHOR = """      fit$y_original <- y_clean
      fit$is_binary <- is_binary_outcome"""

SOFR_SETTINGS_NEW = """      fit$y_original <- y_clean
      fit$is_binary <- is_binary_outcome

      # MERGED APP: what this fit actually used, for the code export.
      fit$fck_settings <- list(
        response     = input$sofr_response,
        predictors   = preds,
        formula      = formula_str,
        family       = pfr_family$family,
        link         = pfr_family$link,
        using_smoothed = !is.null(values$smooth_data),
        n_obs        = length(y_clean))"""

HARMONIC_SETTINGS_ANCHOR = """        include_trend = trend_type != "none",  # For backwards compatibility"""

HARMONIC_SETTINGS_NEW = """        include_trend = trend_type != "none",  # For backwards compatibility
        # MERGED APP: the bounds this run actually used, for the code export.
        fck_settings = list(
          use_bounds = use_bounds, mesor_min = mesor_min, mesor_max = mesor_max,
          amplitude_min = amplitude_min, amplitude_max = amplitude_max,
          A_sat_min = A_sat_min, A_sat_max = A_sat_max,
          tau_min = tau_min, tau_max = tau_max),"""

# C1 continued: three new sections in WaPaa's generate_analysis_code(), so the
# exported script covers the cosinor / FoSR / SoFR fits as well as the fPCA /
# fANOVA / clustering ones it already covered.
#
# The cosinor section emits the app's OWN fit_cosinor() via deparse() rather
# than a hand-written re-implementation.  A re-implementation is a second copy
# of 300 lines of fitting logic that can drift from the one that produced the
# numbers; dumping the real function cannot.
CODE_LIBS_ANCHOR = '''    if(!is.null(values$clustering_results)) {
      add("library(cluster)    # Clustering diagnostics")
      add("library(fda.usc)    # Functional clustering")
    }
    add("")'''

CODE_LIBS_NEW = '''    if(!is.null(values$clustering_results)) {
      add("library(cluster)    # Clustering diagnostics")
      add("library(fda.usc)    # Functional clustering")
    }
    if(!is.null(values$sofr_model)) {
      add("library(refund)     # Scalar-on-function regression")
    }
    if(!is.null(values$harmonic_model) &&
       identical(values$harmonic_model$trend_type, "exp_sat")) {
      add("library(minpack.lm) # Robust non-linear fits (exponential saturation)")
    }
    add("")'''

CODE_SECTIONS_ANCHOR = '''    add("# =============================================================================")
    add("# END OF ANALYSIS CODE")
    add("# =============================================================================")'''

CODE_SECTIONS_NEW = r"""    # ---- SECTION 10: HARMONIC (COSINOR) REGRESSION ----
    # MERGED APP: CIRCAREG had no code generator, so these three sections are
    # new. Settings come from the fit itself (fck_settings), never from the
    # live widgets, so the script matches the model that is on screen.
    if(!is.null(values$harmonic_model)) {
      hm <- values$harmonic_model
      st <- hm$fck_settings
      add("# -----------------------------------------------------------------------------")
      add("# 10. HARMONIC (COSINOR) REGRESSION")
      add("# -----------------------------------------------------------------------------")
      add("")
      add(sprintf("period        <- %s", hm$period))
      add(sprintf("n_harmonics   <- %s", hm$n_harmonics))
      add(sprintf('trend_type    <- "%s"', hm$trend_type))
      add(sprintf("time_vec      <- c(%s)",
                  paste(signif(as.numeric(hm$time_vec), 10), collapse = ", ")))
      add(sprintf("cosinor_input <- %s   # %s",
                  if(isTRUE(hm$using_smoothed)) "smooth_curves" else "raw_data",
                  if(isTRUE(hm$using_smoothed)) "the smoothed curves"
                  else "the raw data (no smoothing had been applied)"))
      add("")
      if(!is.null(st)) {
        add(sprintf("use_bounds    <- %s", isTRUE(st$use_bounds)))
        for(nm in c("mesor_min", "mesor_max", "amplitude_min", "amplitude_max",
                    "A_sat_min", "A_sat_max", "tau_min", "tau_max")) {
          v <- st[[nm]]
          add(sprintf("%-13s <- %s", nm,
                      if(is.null(v) || !is.finite(v)) "NA" else format(v)))
        }
        add("")
      }
      add("# The app's own fitting function, emitted verbatim (not a")
      add("# re-implementation) so this script reproduces its numbers exactly:")
      add(paste("fit_cosinor <-", paste(deparse(fit_cosinor), collapse = "\n")))
      add("")
      if(identical(hm$trend_type, "exp_sat") && exists("fit_cosinor_nonlinear")) {
        add(paste("fit_cosinor_nonlinear <-",
                  paste(deparse(fit_cosinor_nonlinear), collapse = "\n")))
        add("")
      }
      add("# The same per-subject call the app makes:")
      add("cosinor_fits <- lapply(seq_len(nrow(cosinor_input)), function(i) {")
      add("  fit_cosinor(time_vec, cosinor_input[i, ], period = period,")
      add("              n_harmonics = n_harmonics, trend_type = trend_type,")
      if(!is.null(st)) {
        add("              use_bounds = use_bounds,")
        add("              mesor_min = mesor_min, mesor_max = mesor_max,")
        add("              amplitude_min = amplitude_min, amplitude_max = amplitude_max,")
        add("              A_sat_min = A_sat_min, A_sat_max = A_sat_max,")
        add("              tau_min = tau_min, tau_max = tau_max)")
      } else {
        add("              use_bounds = FALSE)")
      }
      add("})")
      add("")
      add("# MESOR / amplitude / acrophase per subject (successful fits only):")
      add("ok <- vapply(cosinor_fits, function(f) isTRUE(f$success), logical(1))")
      add("cosinor_params <- do.call(rbind, lapply(which(ok), function(i) {")
      add("  f <- cosinor_fits[[i]]")
      add("  data.frame(subject = i, mesor = f$mesor,")
      add("             amplitude = f$amplitude[1], acrophase = f$acrophase[1],")
      add("             r_squared = f$r_squared)")
      add("}))")
      add("print(head(cosinor_params))")
      if(!is.null(hm$individual_fits)) {
        add(sprintf("# The app fitted %d of %d subjects successfully.",
                    sum(vapply(hm$individual_fits,
                               function(f) isTRUE(f$success), logical(1))),
                    length(hm$individual_fits)))
      }
      add("")
      add("# NOTE: acrophase is CIRCULAR. Do not average it arithmetically --")
      add("# average the (cos, sin) components, as the app's polar plot does.")
      add("")
    }

    # ---- SECTION 11: FUNCTION-ON-SCALAR REGRESSION ----
    if(!is.null(values$reg_model) && !is.null(values$reg_model$fck_settings)) {
      st <- values$reg_model$fck_settings
      add("# -----------------------------------------------------------------------------")
      add("# 11. FUNCTION-ON-SCALAR REGRESSION (FoSR)")
      add("# -----------------------------------------------------------------------------")
      add("")
      add(sprintf("fosr_predictors <- c(%s)",
                  paste0('"', st$predictors, '"', collapse = ", ")))
      add(sprintf("Y <- %s   # %d subjects x %d time points",
                  if(isTRUE(st$using_smoothed)) "smooth_curves" else "raw_data",
                  st$n_subjects, st$n_time))
      add("keep   <- complete.cases(covariates[, fosr_predictors, drop = FALSE])")
      add("df_reg <- covariates[keep, , drop = FALSE]")
      add("Y      <- Y[keep, , drop = FALSE]")
      add("")
      if(identical(st$method, "OLS_nosmooth")) {
        add("# Pointwise OLS: one regression per time point, in closed form.")
        add("f           <- as.formula(paste('~', paste(fosr_predictors, collapse = ' + ')))")
        add("X           <- model.matrix(f, data = df_reg)")
        add("xtx_inv     <- solve(crossprod(X))")
        add("projector   <- xtx_inv %*% t(X)")
        add("beta_hat    <- projector %*% Y          # the coefficient curves")
        add("fitted_vals <- X %*% beta_hat")
        add("residuals   <- Y - fitted_vals")
        add("n <- nrow(Y); p <- ncol(X)")
        add("sigma2 <- colSums(residuals^2) / (n - p)")
        add("")
        if(isTRUE(st$use_bootstrap)) {
          add(sprintf("# Residual bootstrap for the coefficient bands (B = %d).", st$n_boot))
          add("# The app does not fix a seed, so its bands differ slightly run to run;")
          add("# this script fixes one so the script itself is reproducible.")
          add("set.seed(1)")
          add(sprintf("B <- %d", st$n_boot))
          add("boot_betas <- array(NA, dim = c(B, nrow(beta_hat), ncol(beta_hat)))")
          add("for (b in 1:B) {")
          add("  resid_idx <- sample(seq_len(n), n, replace = TRUE)")
          add("  boot_betas[b, , ] <- projector %*% (fitted_vals + residuals[resid_idx, ])")
          add("}")
          add("boot_ci_lower <- apply(boot_betas, c(2, 3), quantile, 0.025, na.rm = TRUE)")
          add("boot_ci_upper <- apply(boot_betas, c(2, 3), quantile, 0.975, na.rm = TRUE)")
          add("")
        }
      } else {
        add("# Smoothed OLS: one GAM over the (subject, time) long form, with a")
        add("# penalised spline in time and a by-smooth per predictor.")
        add("long_data <- data.frame(")
        add("  value   = as.vector(t(Y)),")
        add("  time    = rep(seq(0, 1, length.out = ncol(Y)), times = nrow(Y)),")
        add("  subject = rep(seq_len(nrow(Y)), each = ncol(Y)))")
        add("for (v in fosr_predictors) long_data[[v]] <- rep(df_reg[[v]], each = ncol(Y))")
        add("gam_formula <- as.formula(paste('value ~ s(time) +',")
        add("  paste(sprintf('s(time, by = %s)', fosr_predictors), collapse = ' + '),")
        add("  '+', paste(fosr_predictors, collapse = ' + ')))")
        add("gam_fit <- mgcv::gam(gam_formula, data = long_data, method = 'REML')")
        add("summary(gam_fit)")
        add("")
        add("# NOTE: the app turns this fit into coefficient curves by prediction")
        add("# contrasts (predict at x = 1 minus predict at x = 0), which is exact")
        add("# for numeric predictors and an approximation for factors.")
        add("")
      }
      add(sprintf("# Method as fitted: %s", values$reg_model$method))
      add("")
    }

    # ---- SECTION 12: SCALAR-ON-FUNCTION REGRESSION ----
    if(!is.null(values$sofr_model) && !is.null(values$sofr_model$fck_settings)) {
      st <- values$sofr_model$fck_settings
      add("# -----------------------------------------------------------------------------")
      add("# 12. SCALAR-ON-FUNCTION REGRESSION (SoFR)")
      add("# -----------------------------------------------------------------------------")
      add("")
      add(sprintf("X_func <- %s   # the curves, as the functional predictor",
                  if(isTRUE(st$using_smoothed)) "smooth_curves" else "raw_data"))
      add(sprintf('y      <- covariates[["%s"]]', st$response))
      if(length(st$predictors)) {
        add(sprintf("sofr_scalars <- c(%s)",
                    paste0('"', st$predictors, '"', collapse = ", ")))
      }
      add("")
      add("keep   <- complete.cases(y) & complete.cases(X_func)")
      add("X_func <- X_func[keep, , drop = FALSE]; y <- y[keep]")
      add("")
      add("pfr_data <- list(X_func = X_func, y = y)")
      if(length(st$predictors)) {
        add("for (v in sofr_scalars) pfr_data[[v]] <- covariates[keep, v]")
      }
      add(sprintf("sofr_fit <- refund::pfr(%s,", st$formula))
      add(sprintf("                        data = pfr_data, family = %s(link = '%s'))",
                  st$family, st$link))
      add("summary(sofr_fit)")
      add("")
      add(sprintf("# Fitted on %d observations; %s family, %s link.",
                  st$n_obs, st$family, st$link))
      add("")
    }

    add("# =============================================================================")
    add("# END OF ANALYSIS CODE")
    add("# =============================================================================")"""


CV_NBASIS_ANCHOR = """          nb <- min(20, n_time - 2)"""

CV_NBASIS_NEW = """          # MERGED APP: follow the user's basis count (CIRCAREG's behaviour) so
          # the recommended smoothing factor applies to the smoothing they will
          # actually run. WaPaa fixed this at 20 regardless.
          nb_user <- if(!is.null(input$smooth_method) && input$smooth_method == "manual")
            input$n_basis_manual else input$n_basis
          if(is.null(nb_user) || !is.finite(nb_user)) nb_user <- 20
          nb <- max(4, min(as.integer(nb_user), n_time - 2))"""

FPCA_FD_ANCHOR = """      if(is.null(values$fd_obj)) {
        n_basis <- min(20, n - 2, max(4, floor(n/3)))
        basis <- create.bspline.basis(c(0, 1), nbasis = n_basis)
        
        # Check if data is already smoothed to avoid double smoothing
        if(!is.null(values$smooth_data)) {
          # Data already smoothed in preprocessing - just create fd representation
          cat("Creating fd_obj from already smoothed data (lambda=0)\\n")
          values$fd_obj <- smooth.basis(time_grid, t(data_to_use), fdPar(basis, 2, 0))$fd
        } else {
          # Raw data - apply default smoothing
          cat("Creating fd_obj from raw data with default smoothing\\n")
          values$fd_obj <- smooth.basis(time_grid, t(data_to_use), basis)$fd
        }
        cat("Created fd_obj with", n_basis, "basis functions\\n")
      }"""

FPCA_FD_NEW = """      if(is.null(values$fd_obj)) {
        # MERGED APP: one rule for this, in FCK/server/04_helpers_fd.R. No
        # smoothing is invented here; an interpolating basis is used and the
        # user is told. (WaPaa silently smoothed onto min(20, n-2) bases.)
        if(!fck_ensure_fd_obj(values)) return()
        cat("Created fd_obj with an interpolating basis (no smoothing)\\n")
      }"""

FANOVA_FD_ANCHOR = """        # CRITICAL: Check if data has already been smoothed in Data Preprocessing
        if(is.null(values$fd_obj)) {
          # Determine which data to use
          data_for_fd <- if(!is.null(values$smooth_data)) {
            cat("Using already smoothed data (no additional smoothing)\\n")
            values$smooth_data
          } else {
            cat("Using raw data (will create fd object)\\n")
            values$data
          }
          
          n_time <- ncol(data_for_fd)
          time_points <- seq(0, 1, length.out = n_time)
          basis <- create.bspline.basis(rangeval = c(0, 1), nbasis = min(20, n_time-2))
          
          # Create fd object WITHOUT additional smoothing (lambda = 0 for already smoothed data)
          if(!is.null(values$smooth_data)) {
            # Data already smoothed - just create fd representation with no penalty
            values$fd_obj <- smooth.basis(time_points, t(data_for_fd), fdPar(basis, 2, 0))$fd
          } else {
            # Raw data - apply default smoothing
            values$fd_obj <- smooth.basis(time_points, t(data_for_fd), basis)$fd
          }
        }"""

FANOVA_FD_NEW = """        # CRITICAL: Check if data has already been smoothed in Data Preprocessing.
        # MERGED APP: when it has not, fck_ensure_fd_obj() builds an
        # INTERPOLATING representation and says so, instead of quietly
        # smoothing onto min(20, n_time - 2) basis functions as WaPaa did.
        # Nothing here re-smooths data that the preprocessing step already
        # smoothed: values$fd_obj is reused as-is.
        if(is.null(values$fd_obj)) {
          if(!fck_ensure_fd_obj(values)) return()
        }"""

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
     [("W", 1406, 1494, "smoothing fit statistics printout", "time_axis"),
      ("C", 913, 930, "compact smoothing fit-metrics panel", None),
      ("W", 2862, 3016, "interactive smoothed-curve plot + curve selection", None)]),

    # ---- server: shared smoothing diagnostics ------------------------------
    # Both apps had this section; WaPaa's is a strict superset of CIRCAREG's
    # (it adds stratified CV folds and the GCV-vs-n-basis sweep) and its
    # reactiveValues names are the ones the merged app uses, so WaPaa's is
    # carried across whole and CIRCAREG's duplicate is dropped.
    ("server/30_diagnostics.R", None,
     [("W", 2231, 2860, "GAM REML, REML profile, CV, n-basis sweep", "cv_nbasis")]),

    # ---- server: analysis families ----------------------------------------
    ("server/40_fpca.R", None,
     [("W", 3018, 4423, "group summary, fPCA/warping analysis + outputs", "fpca_fd"),
      ("W", 6229, 6522, "warping / alignment / landmark plots", "landmark")]),
    ("server/50_fanova.R", None,
     [("W", 4425, 6227, "group UIs, fANOVA (between + repeated measures)", "fanova_fd"),
      ("W", 6524, 6976, "post-hoc pairwise outputs", None)]),
    ("server/60_clustering.R", None,
     [("W", 6978, 8912, "functional clustering, optimisation, DCF, outputs", None)]),
    ("server/70_fosr.R", None,
     [("C", 1658, 2242, "Function-on-Scalar regression", "fosr_settings")]),
    ("server/71_sofr.R", None,
     [("C", 2244, 2877, "Scalar-on-Function regression", "sofr_settings")]),
    ("server/72_harmonic.R", None,
     [("C", 2879, 7190, "cosinor core, harmonic regression + outputs", "harmonic")]),
    ("server/73_cosinor_pairwise.R", None,
     [("C", 7264, 7859, "cosinor pairwise group tests; ids prefixed hp_", "hp")]),
    ("server/90_export.R", None,
     [("W", 8914, 9929, "all WaPaa exports + the reproducible-code generator", "codegen"),
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
                chunk = patch(chunk, HARMONIC_SETTINGS_ANCHOR,
                              HARMONIC_SETTINGS_NEW, rel)
            elif transform == "codegen":
                chunk = patch(chunk, CODE_LIBS_ANCHOR, CODE_LIBS_NEW, rel)
                chunk = patch(chunk, CODE_SECTIONS_ANCHOR, CODE_SECTIONS_NEW, rel)
            elif transform == "fosr_settings":
                chunk = patch(chunk, FOSR_SETTINGS_ANCHOR, FOSR_SETTINGS_NEW, rel)
            elif transform == "sofr_settings":
                chunk = patch(chunk, SOFR_SETTINGS_ANCHOR, SOFR_SETTINGS_NEW, rel)
            elif transform == "cv_nbasis":
                chunk = patch(chunk, CV_NBASIS_ANCHOR, CV_NBASIS_NEW, rel)
            elif transform == "fpca_fd":
                chunk = patch(chunk, FPCA_FD_ANCHOR, FPCA_FD_NEW, rel)
            elif transform == "fanova_fd":
                chunk = patch(chunk, FANOVA_FD_ANCHOR, FANOVA_FD_NEW, rel)
            elif transform == "time_axis":
                chunk = patch(chunk, TIME_AXIS_ANCHOR, TIME_AXIS_NEW, rel)
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
