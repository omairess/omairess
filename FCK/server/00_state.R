# ==============================================================================
# server/00_state.R — THE shared state bus
#
# Union of the two source apps' reactiveValues (WaPaa1_3.R lines 1125-1147 and
# CIRCAREG.R lines 684-700).  Where both apps stored the same thing under
# different names, the merged app keeps ONE name and tools/port_fck.py renames
# the ported code accordingly:
#
#   CIRCAREG                     merged app
#   ---------------------------  ---------------------------
#   gam_reml_result              gam_reml_fit
#   reml_profile_result          reml_profile
#   cv_result                    cv_results
#   pairwise_results (cosinor)   hp_pairwise_results
#
# The first three are smoothing diagnostics: the merged app runs WaPaa's
# version of that section (a strict superset), so CIRCAREG's names disappear
# with its duplicate code.  The fourth moves because WaPaa's post-hoc tests
# already own `pairwise_results`.
#
# THE DATA CONTRACT every analysis tab reads:
#   values$data         numeric matrix, subjects x time points, as selected
#   values$smooth_data  the same matrix after the shared smoothing step
#   values$fd_obj       fda object for the smoothed curves, on a 0-1 range
#                       (fPCA / warping / fANOVA / clustering assume 0-1)
#   values$time_labels  original column names, in file order
#   values$time_numeric column indices 1:n_time (WaPaa's plotting x axis)
#   values$time_clock   real clock hours from the column names, or NULL
#   values$covariates   scalar variables, ORIGINAL types (predictors/response)
#   values$group_variables / $selected_group_vars / $group_labels
#                       the same scalar variables as factors, for grouping
# ==============================================================================

values <- reactiveValues(
  # ---- shared: data ---------------------------------------------------------
  raw_df                = NULL,   # raw imported data frame
  uploaded_data         = NULL,   # same, kept for the RM-ANOVA variable pickers
  data                  = NULL,   # numeric analysis matrix (subjects x time)
  time_labels           = NULL,   # original column names for the time points
  time_numeric          = NULL,   # WaPaa's plotting x coordinates: 1:n_time.
                                  #   NOT clock times — see 03_helpers_clock.R
  time_clock            = NULL,   # real clock hours parsed from the column
                                  #   names, or NULL when they cannot be
                                  #   trusted. This is the one to use when an
                                  #   analysis needs actual elapsed time.
  covariates            = NULL,   # scalar variables, original types
  group_variables       = NULL,   # the same variables coerced to factors
  selected_group_vars   = NULL,   # their names
  group_labels          = NULL,   # primary grouping factor (first selected)

  # ---- shared: smoothing ----------------------------------------------------
  smooth_data           = NULL,
  fd_obj                = NULL,
  smooth_fit_metrics    = NULL,   # R2 / RMSE / EDF / GCV per subject
  smoothing_avg_metrics = NULL,   # averaged metrics, WaPaa display format

  # ---- shared: smoothing diagnostics ---------------------------------------
  gam_reml_fit          = NULL,
  gam_data_label        = NULL,
  reml_profile          = NULL,
  cv_results            = NULL,
  diagnostic_lambda     = NULL,
  nbasis_gcv            = NULL,

  # ---- F: functional PCA / warping / ANOVA / clustering (from WaPaa) --------
  pca_results           = NULL,
  warping_results       = NULL,
  landmarks             = list(),
  landmark_points       = data.frame(x = numeric(), y = numeric()),
  fanova_results        = NULL,
  fanova_selected_groups = NULL,
  pairwise_results      = NULL,   # post-hoc tests after functional ANOVA
  clustering_results    = NULL,
  cluster_optimization  = NULL,
  selected_curve        = NULL,

  # ---- C: functional / circadian regression (from CIRCAREG) ----------------
  reg_model             = NULL,   # function-on-scalar
  sofr_model            = NULL,   # scalar-on-function
  harmonic_model        = NULL,   # cosinor
  hp_pairwise_results   = NULL,   # pairwise tests on cosinor parameters
  hp_pairwise_param     = NULL,
  hp_pairwise_correction = NULL
)

# Everything downstream of the shared data step is invalidated whenever the
# data step re-runs.  Both source apps did this partially and for their own
# analyses only; here it is one list so a new selection cannot leave a stale
# cosinor fit sitting next to a fresh fPCA.
fck_reset_analyses <- function(values, keep_smoothing = FALSE) {
  if (!keep_smoothing) {
    values$smooth_data <- NULL
    values$fd_obj <- NULL
    values$smooth_fit_metrics <- NULL
    values$smoothing_avg_metrics <- NULL
  }
  values$pca_results <- NULL
  values$warping_results <- NULL
  values$fanova_results <- NULL
  values$pairwise_results <- NULL
  values$clustering_results <- NULL
  values$cluster_optimization <- NULL
  values$reg_model <- NULL
  values$sofr_model <- NULL
  values$harmonic_model <- NULL
  values$hp_pairwise_results <- NULL
  invisible(NULL)
}
