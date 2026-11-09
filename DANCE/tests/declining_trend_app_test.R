# ==============================================================================
# tests/declining_trend_app_test.R -- a saturating trend that FALLS
#
# A participant whose response falls over the recording (a negative linear
# slope) used to get a start value that carried lm()'s coefficient name, so
# nls called the parameter "A_sat.t_shifted". With tau held the fit was thrown
# away ("Result extraction failed"), with tau estimated A_sat came back NA --
# either way the participant lost the saturating trend, and the warning said
# "need 8+ valid points" for participants with every point. Here half the
# participants fall: every one must be fitted, with a negative A, under both
# tau settings, through the app and the exported script; a flat participant is
# reported as having no variation, not as lacking points.
# Run with:   Rscript tests/declining_trend_app_test.R      (from the DANCE directory)
# ==============================================================================
src <- readLines(if (dir.exists("tests")) "tests/structures_app_test.R" else "DANCE/tests/structures_app_test.R")
eval(parse(text = src[1:(grep("^# 24 participants", src) - 1)]))
set.seed(52)
clock <- c(8, 10, 12, 14, 16, 18, 20, 22, 0, 2, 4, 6, 8)
tt <- seq(0, 24, by = 2)
n <- 24; grp <- rep(c("0", "1"), each = n / 2)
sgn <- rep(c(1, -1), n / 2)                       # every other participant FALLS
Y <- t(vapply(seq_len(n), function(i)
  40 + rnorm(1, 0, 4) + sgn[i] * (20 + rnorm(1, 0, 3)) * (1 - exp(-tt / 10)) +
    8 * cos(2 * pi * (tt + 8 - 4) / 24) + rnorm(length(tt), 0, 3), numeric(length(tt))))
Y[7, ] <- 15                                       # one participant with no variation at all
d <- data.frame(ID = sprintf("P%02d", seq_len(n)), sex = grp, round(Y, 2))
names(d)[-(1:2)] <- c(paste0("DZ_", clock[1:12]), "DZ_8_SD")
csv <- tempfile(fileext = ".csv"); write.csv(d, csv, row.names = FALSE, na = "")
server <- function(input, output, session) for (f in server_files) dance_source(f, envir = environment())

testServer(server, {
  session$setInputs(datafile = list(name = "dz.csv", datapath = csv), header = TRUE, sep = ",",
                    data_format = "wide", load_data = 1)
  session$setInputs(sel_data_vars = names(values$raw_df)[-(1:2)], sel_cov_vars = "sex", subject_id_var = "ID",
                    apply_selection = 1)
  k <- 0L
  for (ts_mode in c("fixed", "estimate")) {
    k <- k + 1L
    session$setInputs(harmonic_approach = "two_stage", harmonic_structure = "combined", harmonic_trend_type = "exp_sat",
                      harmonic_period = 24, n_harmonics = 2, harmonic_time_var = "_shared_",
                      harmonic_group_var = "sex", harmonic_design = "between", harmonic_between_var = "sex",
                      harmonic_ts_tau_mode = ts_mode, harmonic_traj_tau = 10, harmonic_tau_lo = 1, harmonic_tau_hi = 72,
                      harmonic_origin = "first", run_harmonic = k)
    mod <- values$harmonic_model
    ok <- vapply(mod$individual_fits, function(f) isTRUE(f$success), logical(1))
    fall <- which(sgn < 0 & seq_len(n) != 7L)
    chk(all(ok[-7]) && !ok[7],
        sprintf("tau %s: every participant with variation is fitted, falling ones included", ts_mode),
        sprintf("tau %s: not fitted: %s", ts_mode, paste(which(!ok), collapse = ", ")))
    A <- vapply(mod$individual_fits[fall], function(f) f$trend_params$A_sat$coef %||% NA_real_, numeric(1))
    chk(all(is.finite(A)) && mean(A < 0) >= 0.8,
        sprintf("tau %s: the falling participants have A estimated, %d of %d negative", ts_mode, sum(A < 0), length(A)),
        sprintf("tau %s: A of the falling participants: %s", ts_mode, paste(round(A, 2), collapse = ", ")))
    msg <- dance_ts_failure_text(Filter(Negate(is.null), lapply(seq_len(n), function(i)
      if (!ok[i]) list(subject = i, n_valid = sum(is.finite(Y[i, ]))))), mod$Y, 8L, mixed = FALSE)
    chk(grepl("no variation, every value equal (nothing to fit): S7", msg, fixed = TRUE) && !grepl("too few", msg),
        sprintf("tau %s: the flat participant is reported as having no variation", ts_mode),
        sprintf("tau %s: the failure message is wrong: %s", ts_mode, msg))
    # the exported script fits the same participants and reproduces the outcomes
    code <- generate_analysis_code()
    env <- new.env(parent = globalenv()); env$dance_input_data <- values$data
    ran <- tryCatch({ suppressWarnings(suppressMessages(capture.output(eval(parse(text = code), envir = env)))); TRUE },
                    error = function(e) { cat("   script error:", conditionMessage(e), "\n"); FALSE })
    so <- env$ts_outcomes; oc <- harmonic_outcomes()
    chk(ran && !is.null(so) && isTRUE(all.equal(so$estimates$estimate, oc$estimates$estimate, tolerance = 1e-6)) &&
          sum(vapply(env$ts_fit$fits, function(f) isTRUE(f$success), logical(1))) == sum(ok),
        sprintf("tau %s: the exported script fits the same participants and reproduces every estimate", ts_mode),
        sprintf("tau %s: the script differs from the app", ts_mode))
  }
})
cat(sprintf("\n%s\n", if (failures == 0) "ALL CHECKS PASSED" else sprintf("%d CHECK(S) FAILED", failures)))
if (failures > 0) quit(status = 1)
