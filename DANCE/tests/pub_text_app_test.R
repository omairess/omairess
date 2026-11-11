# Publication-ready text in the analysis report (server/08p_helpers_pubtext.R):
# the section exists for both approaches; every number in it is one of the
# outcome table's; pairwise differences are stated positive and oriented; the
# wording has no doubled nouns and every sentence starts with a capital.
src <- readLines("tests/structures_app_test.R")
eval(parse(text = src[1:(grep("^testServer", src) - 1)]))
fails <- 0L
chk <- function(ok, what) { cat(if (isTRUE(ok)) "PASS" else { fails <<- fails + 1L; "FAIL" }, what, "\n") }
for (AP in c("mixed", "two_stage")) testServer(server, {
  session$setInputs(datafile = list(name = "kss.csv", datapath = csv), header = TRUE, sep = ",", data_format = "wide", load_data = 1)
  session$setInputs(sel_data_vars = names(values$raw_df)[-(1:2)], sel_cov_vars = "group", subject_id_var = "ID", apply_selection = 1)
  session$setInputs(harmonic_approach = AP, harmonic_structure = "combined", harmonic_trend_type = "exp_sat",
                    harmonic_period = 24, n_harmonics = 1, harmonic_time_var = "_shared_", harmonic_group_var = "group",
                    harmonic_design = "between", harmonic_between_var = "group", harmonic_df_method = "kr",
                    harmonic_traj_adjust = "none", harmonic_traj_resid = "independent", harmonic_traj_tau = 14,
                    harmonic_traj_tau_mode = "fixed", harmonic_ts_tau_mode = "fixed", harmonic_origin = "first", run_harmonic = 1)
  md <- dance_apa_report(values, input)
  i <- grep("^## Publication-ready text", md); j <- grep("^## Reproducibility", md)
  chk(length(i) == 1L && length(j) == 1L && i < j, paste(AP, "section present before Reproducibility"))
  pub <- md[i:(j - 1)]; txt <- paste(pub, collapse = "\n")
  for (hd in c("### Statistical Analysis", "### Results", "#### Overall response level",
               "#### Non-periodic temporal change", "#### Rhythmic component", "#### Summary of findings", "#### Table 1"))
    chk(any(pub == hd), paste(AP, "heading", hd))
  o <- harmonic_outcomes(); pw <- harmonic_pairwise_all()
  f2 <- function(x) formatC(x, format = "f", digits = 2)
  # group estimates with their CIs, as the outcome table has them
  e <- o$estimates; e <- e[e$outcome == "fitted_mean", ]
  chk(nrow(e) == 2L, paste(AP, "two level estimates in the outcome table"))
  for (k in seq_len(nrow(e)))
    chk(grepl(sprintf("(%s, 95%% CI [%s, %s])", f2(e$estimate[k]), f2(e$lo[k]), f2(e$hi[k])), txt, fixed = TRUE),
        paste(AP, "level estimate and CI of", e$cell[k]))
  # the significant level contrast: positive, the larger group named first
  r <- pw[pw$quantity == "fitted_mean", ][1, ]
  hi_first <- r$estimate > 0
  ab <- strsplit(r$comparison, " vs ")[[1]]; big <- if (hi_first) ab[1] else ab[2]; small <- setdiff(ab, big)
  chk(grepl(sprintf("the %s group had a higher time-averaged level than the %s group (difference = %s", big, small,
                    f2(abs(r$estimate))), txt, fixed = TRUE), paste(AP, "level contrast oriented and positive"))
  chk(!grepl("difference = -", txt, fixed = TRUE), paste(AP, "no negative stated difference"))
  chk(!grepl("group groups", txt, fixed = TRUE), paste(AP, "no 'group groups'"))
  body <- pub[!grepl("^(#|\\||\\*|$)", pub)]
  sent <- unlist(strsplit(body, "(?<=\\.) ", perl = TRUE))
  chk(!any(grepl("^[a-z]", trimws(sent))), paste(AP, "every sentence starts with a capital"))
  chk(!any(grepl("[ \t]$", pub)), paste(AP, "no trailing whitespace"))
  chk(grepl("*F*(", paste(pub[grep("joint test of the harmonic", pub)], collapse = ""), fixed = TRUE),
      paste(AP, "rhythm joint test stated with its F"))
  chk(grepl("^\\| Parameter \\| ctrl \\(\\*n\\* = 12\\) \\| sd \\(\\*n\\* = 12\\) \\|$", pub[grep("^\\| Parameter", pub)]),
      paste(AP, "Table 1 header uses the group labels"))
})
# numeric levels and three groups (as in a factor coded 1/2/3): "Group 1",
# "the age groups", lists joined with "or"
d3 <- read.csv(csv); d3$age_group <- rep(1:3, length.out = nrow(d3)); d3$group <- NULL
csv3 <- tempfile(fileext = ".csv")
write.csv(d3[, c("ID", "age_group", setdiff(names(d3), c("ID", "age_group")))], csv3, row.names = FALSE, na = "")
testServer(server, {
  session$setInputs(datafile = list(name = "kss.csv", datapath = csv3), header = TRUE, sep = ",", data_format = "wide", load_data = 1)
  session$setInputs(sel_data_vars = names(values$raw_df)[-(1:2)], sel_cov_vars = "age_group", subject_id_var = "ID", apply_selection = 1)
  session$setInputs(harmonic_approach = "mixed", harmonic_structure = "combined", harmonic_trend_type = "exp_sat",
                    harmonic_period = 24, n_harmonics = 1, harmonic_time_var = "_shared_", harmonic_group_var = "age_group",
                    harmonic_design = "between", harmonic_between_var = "age_group", harmonic_df_method = "kr",
                    harmonic_traj_adjust = "none", harmonic_traj_resid = "independent", harmonic_traj_tau = 14,
                    harmonic_traj_tau_mode = "fixed", harmonic_origin = "first", run_harmonic = 1)
  md <- dance_apa_report(values, input); txt <- paste(md, collapse = "\n")
  chk(grepl("| Parameter | Group 1 (*n* = 8) | Group 2 (*n* = 8) | Group 3 (*n* = 8) |", txt, fixed = TRUE), "3 groups: numeric labels in Table 1")
  chk(grepl("between the age groups", txt, fixed = TRUE) && !grepl("between age group,", txt, fixed = TRUE), "3 groups: 'the age groups'")
  chk(grepl("Groups 1 and 3 \\(\\*p\\* = [.0-9]+\\), or Groups 2 and 3", txt), "3 groups: pairwise list joined with 'or'")
  chk(grepl("*Model-based parameter estimates by age group.*", txt, fixed = TRUE), "3 groups: Table 1 title")
})
cat(if (fails) sprintf("%d FAILED\n", fails) else "ALL PASS\n"); quit(status = fails > 0L)
