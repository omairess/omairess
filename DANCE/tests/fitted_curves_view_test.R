# ==============================================================================
# tests/fitted_curves_view_test.R -- the Fitted curves and trajectory plots
#
#   * one participant's raw data belong to the participant's legend entry, so
#     hiding their group's curve does not hide them;
#   * raw points are drawn large enough to see (size 8);
#   * the components (each harmonic and the trend) can be switched on under the
#     mixed-effects approach too, each its own legend entry;
#   * the bands can be switched off (Estimated trajectories: "None").
# Run with:   Rscript tests/fitted_curves_view_test.R      (from the DANCE directory)
# ==============================================================================
src <- readLines(if (dir.exists("tests")) "tests/structures_app_test.R" else "DANCE/tests/structures_app_test.R")
eval(parse(text = src[1:(grep("^testServer", src) - 1)]))
built <- function(p) plotly::plotly_build(p)$x$data

testServer(server, {
  session$setInputs(datafile = list(name = "kss.csv", datapath = csv), header = TRUE, sep = ",",
                    data_format = "wide", load_data = 1)
  session$setInputs(sel_data_vars = names(values$raw_df)[-(1:2)], sel_cov_vars = "group", subject_id_var = "ID",
                    apply_selection = 1)
  session$setInputs(harmonic_approach = "mixed", harmonic_structure = "combined", harmonic_trend_type = "exp_sat",
                    harmonic_period = 24, n_harmonics = 2, harmonic_time_var = "_shared_",
                    harmonic_group_var = "group", harmonic_design = "between", harmonic_between_var = "group",
                    harmonic_df_method = "satterthwaite", harmonic_traj_resid = "independent",
                    harmonic_traj_tau = 14, harmonic_traj_tau_mode = "fixed", harmonic_origin = "first",
                    run_harmonic = 1)
  # ONE participant with raw data
  session$setInputs(harmonic_subject_select = "P03", harmonic_show_data = TRUE, harmonic_raw_style = "points",
                    harmonic_show_ci = TRUE, harmonic_show_components = FALSE)
  d <- built(harmonic_fit_plot_build())
  raw <- Filter(function(t) identical(t$mode, "markers"), d)
  chk(length(raw) == 1 && identical(raw[[1]]$legendgroup, "participant") && isTRUE(raw[[1]]$showlegend),
      "one participant's data are their own legend entry (hiding the group keeps them)",
      sprintf("the participant's data are in legend group '%s'", if (length(raw)) raw[[1]]$legendgroup else "none"))
  chk(length(raw) == 1 && isTRUE(raw[[1]]$marker$size >= 8), "raw points are size 8",
      sprintf("raw point size %s", if (length(raw)) raw[[1]]$marker$size else NA))
  # COMPONENTS under the mixed approach
  session$setInputs(harmonic_subject_select = "mean", harmonic_show_data = FALSE, harmonic_show_components = TRUE)
  d <- built(harmonic_fit_plot_build())
  nm <- vapply(d, function(t) as.character(t$name %||% ""), character(1))
  chk(sum(grepl(": H1 (24", nm, fixed = TRUE)) == 2 && sum(grepl(": H2 (12", nm, fixed = TRUE)) == 2 &&
        sum(grepl("saturat", nm, ignore.case = TRUE)) == 2,
      "mixed: each group's H1, H2 and trend components are drawn, each its own legend entry",
      sprintf("mixed components: %s", paste(nm, collapse = " | ")))
  session$setInputs(harmonic_show_components = FALSE, harmonic_show_ci = FALSE)
  d <- built(harmonic_fit_plot_build())
  chk(!any(vapply(d, function(t) identical(t$fill, "toself") || identical(t$fill, "tonexty"), logical(1))),
      "Fitted curves: the bands switch off", "Fitted curves: a band is still drawn")
  # ESTIMATED TRAJECTORIES: band "none"
  session$setInputs(harmonic_traj_component = "full", harmonic_traj_band = "none", harmonic_traj_annot = character(0),
                    harmonic_traj_raw = FALSE)
  d <- built(harmonic_traj_curves_build())
  chk(!any(vapply(d, function(t) identical(t$fill, "toself") || identical(t$fill, "tonexty"), logical(1))) &&
        sum(vapply(d, function(t) identical(t$mode, "lines"), logical(1))) >= 2,
      "Estimated trajectories: 'None' draws the curves without bands", "Estimated trajectories: bands still drawn")
  session$setInputs(harmonic_traj_band = "pointwise")
  d <- built(harmonic_traj_curves_build())
  chk(any(vapply(d, function(t) identical(t$fill, "toself") || identical(t$fill, "tonexty"), logical(1))),
      "Estimated trajectories: 'Pointwise' draws the bands again", "Estimated trajectories: no band with Pointwise")
})
cat(sprintf("\n%s\n", if (failures == 0) "ALL CHECKS PASSED" else sprintf("%d CHECK(S) FAILED", failures)))
if (failures > 0) quit(status = 1)
