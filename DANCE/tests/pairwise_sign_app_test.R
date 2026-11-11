# ==============================================================================
# tests/pairwise_sign_app_test.R -- one sign convention for every pairwise row
#
# "a vs b" is a minus b, for every quantity and both approaches: the level and
# trend quantities, the amplitude and the acrophase. The mixed amplitude and
# acrophase rows used to be b minus a (the helpers' own convention), so a
# table could show a positive difference for "1 vs 2" while group 2 had the
# larger amplitude. Each row is checked against the groups' own estimates.
# Run with:   Rscript tests/pairwise_sign_app_test.R      (from the DANCE directory)
# ==============================================================================
src <- readLines(if (dir.exists("tests")) "tests/structures_app_test.R" else "DANCE/tests/structures_app_test.R")
eval(parse(text = src[1:(grep("^testServer", src) - 1)]))

testServer(server, {
  session$setInputs(datafile = list(name = "kss.csv", datapath = csv), header = TRUE, sep = ",",
                    data_format = "wide", load_data = 1)
  session$setInputs(sel_data_vars = names(values$raw_df)[-(1:2)], sel_cov_vars = "group", subject_id_var = "ID",
                    apply_selection = 1)
  k <- 0L
  for (ap in c("mixed", "two_stage")) {
    k <- k + 1L
    session$setInputs(harmonic_approach = ap, harmonic_structure = "combined", harmonic_trend_type = "exp_sat",
                      harmonic_period = 24, n_harmonics = 1, harmonic_time_var = "_shared_",
                      harmonic_group_var = "group", harmonic_design = "between", harmonic_between_var = "group",
                      harmonic_df_method = "satterthwaite", harmonic_traj_resid = "independent",
                      harmonic_traj_tau = 14, harmonic_traj_tau_mode = "fixed", harmonic_ts_tau_mode = "fixed",
                      harmonic_origin = "first", run_harmonic = k)
    oc <- harmonic_outcomes(); es <- oc$estimates
    val <- function(q, cl) es$estimate[es$outcome == q & es$cell == cl][1]
    for (q in c("fitted_mean", "change", "amplitude1", "phase1")) {
      r <- if (ap == "mixed") dance_traj_contrasts(harmonic_traj(), q, adjust = "none", df_method = "satterthwaite")
           else { ts <- harmonic_ts(); dance_ts_pairwise(ts$df, ts$mod, q, adjust = "none") }
      tb <- r$table
      want <- val(q, tb$cell1[1]) - val(q, tb$cell2[1])
      # two-stage amplitude and acrophase rows compare the participants' OWN
      # values (Welch's t on the amplitudes, Watson-Williams on the phases), not
      # the population-mean cosinor's vector-mean amplitude in the group table
      if (ap == "two_stage" && q %in% c("amplitude1", "phase1")) {
        df <- harmonic_ts()$df; g <- as.character(df$group)
        if (q == "amplitude1") want <- mean(df$amplitude_1[g == tb$cell1[1]]) - mean(df$amplitude_1[g == tb$cell2[1]])
        else {
          cm <- function(x) atan2(mean(sin(x)), mean(cos(x)))
          d <- cm(df$acrophase_rad_1[g == tb$cell1[1]]) - cm(df$acrophase_rad_1[g == tb$cell2[1]])
          want <- atan2(sin(d), cos(d)) * 24 / (2 * pi)
        }
      }
      if (q == "phase1") want <- ((want + 12) %% 24) - 12      # the shortest arc
      tol <- if (q %in% c("amplitude1", "phase1") && ap == "mixed") 0.05 * max(1, abs(want)) else 1e-6 * max(1, abs(want))
      chk(isTRUE(r$ok) && isTRUE(abs(tb$estimate[1] - want) < tol) &&
            (is.na(tb$lo[1]) || (tb$lo[1] <= tb$estimate[1] && tb$estimate[1] <= tb$hi[1])),
          sprintf("%s, %s: '%s vs %s' = %.3f = %s minus %s (%.3f)", ap, q, tb$cell1[1], tb$cell2[1],
                  tb$estimate[1], tb$cell1[1], tb$cell2[1], want),
          sprintf("%s, %s: '%s vs %s' = %s but %s minus %s is %.3f", ap, q, tb$cell1[1], tb$cell2[1],
                  format(tb$estimate[1]), tb$cell1[1], tb$cell2[1], want))
    }
  }
})
cat(sprintf("\n%s\n", if (failures == 0) "ALL CHECKS PASSED" else sprintf("%d CHECK(S) FAILED", failures)))
if (failures > 0) quit(status = 1)
