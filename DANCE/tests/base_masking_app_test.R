# ==============================================================================
# tests/base_masking_app_test.R -- a package in the user's session masks a base
# function; the app must still use base R
#
# A package attached BEFORE the app starts can export a function with a base
# name (rsat exports `sub`). Every sub() in the app then reached it, and the
# cosinor summary failed with "subscript out of bounds". app.R now binds base
# R's version of whatever such a package masks, where the server functions
# look first -- except for DANCE's own packages and their dependencies (Matrix,
# attached with lme4, extends t(), solve(), as.matrix() ... and must keep them).
#
# Also: column names ending in _24 (midnight) and a second day marked _SD give
# clock times (10, ..., 22, 24, 2, ..., 8, 10_SD, ..., 20_SD -> 36 h).
#
# Loads app.R the way shiny::runApp does, attaches a masking `sub` first, and
# runs a mixed-effects cosinor on a small two-group sample through the summary.
#
# Run with:   Rscript tests/base_masking_app_test.R      (from the DANCE directory)
# ==============================================================================
.libPaths(c("~/Rlib", .libPaths()))
`%||%` <- function(a, b) if (is.null(a)) b else a
failures <- 0L
chk <- function(cond, good, bad) if (isTRUE(cond)) cat("ok  :", good, "\n") else { cat("FAIL:", bad, "\n"); failures <<- failures + 1L }
app_dir <- if (file.exists("app.R")) "." else "DANCE"

suppressPackageStartupMessages(library(Matrix))   # a DANCE dependency, attached by the user first
mask <- new.env()
assign("sub", function(...) stop("subscript out of bounds (the masking sub was called)"), envir = mask)
attach(mask, name = "package:fakemask", warn.conflicts = FALSE)

# a small two-group sample: 22 participants, 18 two-hourly samples over 36 h
set.seed(3)
hrs <- c(10, 12, 14, 16, 18, 20, 22, 24, 2, 4, 6, 8)
labs <- c(paste0("CORT_", hrs), paste0("CORT_", c(10, 12, 14, 16, 18, 20), "_SD"))
tt <- seq(0, by = 2, length.out = 18); clock <- (10 + tt) %% 24
grp <- rep(0:1, each = 11)
Y <- t(vapply(seq_along(grp), function(i)
  4 + rnorm(1, 0, 1) + 0.6 * grp[i] + (2 + rnorm(1, 0, 0.4)) * cos(2 * pi * (clock - 8) / 24) +
    rnorm(18, 0, 0.8), numeric(18)))
Y[3, 6] <- NA
d <- data.frame(ID = sprintf("K%02d", seq_along(grp)), Group_short_long = grp, round(Y, 2))
names(d)[-(1:2)] <- labs
csv_path <- tempfile(fileext = ".csv"); write.csv(d, csv_path, row.names = FALSE, na = "")

app_env <- new.env(parent = globalenv())
msgs <- character(0)
withCallingHandlers(invisible(sys.source(file.path(app_dir, "app.R"), envir = app_env, chdir = TRUE)),
                    message = function(m) { msgs <<- c(msgs, conditionMessage(m)); invokeRestart("muffleMessage") })
chk(identical(names(app_env$DANCE_BASE_RESTORED), "sub"),
    "app.R binds base sub() over the masking package, and nothing of Matrix's",
    sprintf("restored: %s", paste(names(app_env$DANCE_BASE_RESTORED), collapse = ", ")))
chk(any(grepl("using base R's sub() (masked by fakemask)", msgs, fixed = TRUE)),
    "it says so at startup", "no startup message about the masked function")

shiny::testServer(app_env$server, {
  session$setInputs(datafile = list(name = "cort.csv", datapath = csv_path), header = TRUE, sep = ",",
                    data_format = "wide", load_data = 1)
  session$setInputs(sel_data_vars = labs, sel_cov_vars = "Group_short_long", subject_id_var = "_none_",
                    apply_selection = 1)
  chk(isTRUE(all.equal(dance_cumulative_hours(values$time_labels), tt)),
      "_24 and _SD columns give 36 h of clock time", sprintf("clock: %s", paste(values$time_clock, collapse = ",")))
  session$setInputs(harmonic_approach = "mixed", harmonic_period = 24, n_harmonics = 1,
                    harmonic_trend_type = "none", harmonic_time_var = "_shared_",
                    harmonic_design = "between", harmonic_between_var = "Group_short_long",
                    harmonic_df_method = "kr", harmonic_traj_adjust = "holm",
                    harmonic_traj_what = "intercept", harmonic_traj_resid = "independent",
                    harmonic_tau_lo = 2, harmonic_tau_hi = 48, harmonic_traj_tau = 18,
                    harmonic_traj_tau_mode = "fixed", run_harmonic = 1)
  chk(isTRUE(values$harmonic_model$traj$ok), "the mixed cosinor fits (lme4 with Matrix intact)",
      sprintf("the fit failed: %s", values$harmonic_model$traj$message %||% ""))
  sm <- tryCatch(output$harmonic_summary, error = function(e) paste("ERROR", conditionMessage(e)))
  chk(!grepl("^ERROR", sm) && grepl("Assumptions", sm), "the summary, with its assumptions, prints",
      sprintf("summary: %s", substr(sm, 1, 160)))
  x <- tryCatch(paste(unlist(output$harmonic_traj_primary), collapse = " "), error = function(e) "ERROR")
  chk(grepl("Kenward-Roger", x), "the group test is computed", "the group test failed")
})
detach("package:fakemask")
cat(sprintf("\n%s\n", if (failures == 0) "ALL CHECKS PASSED" else sprintf("%d CHECK(S) FAILED", failures)))
if (failures > 0) quit(status = 1)
