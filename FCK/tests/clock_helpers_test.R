# ==============================================================================
# tests/clock_helpers_test.R — the clock-time parsing behind two features
#
# server/03_helpers_clock.R decides (a) whether the cosinor tab's "use shared
# clock times" option is allowed to fire, and (b) what the opt-in real-time
# smoothing uses as its argument. Both are silently wrong if the parsing is
# wrong — a column counter passed off as hours looks perfectly reasonable on a
# plot — so the shapes that matter are pinned here.
#
# Run with:   Rscript tests/clock_helpers_test.R      (from the FCK directory)
# Needs no packages at all.
# ==============================================================================

app_dir <- if (dir.exists("server")) "." else "FCK"
eval(parse(file.path(app_dir, "server/01_helpers_time.R"), encoding = "UTF-8"))
eval(parse(file.path(app_dir, "server/03_helpers_clock.R"), encoding = "UTF-8"))

failures <- 0
chk <- function(label, labs, expect_hours, expect_uneven, expect_cum = NULL) {
  h <- fck_clock_hours(labs)
  u <- fck_spacing_is_uneven(labs)
  cum <- fck_cumulative_hours(labs)
  bad <- character(0)
  if (!identical(is.null(h), is.null(expect_hours))) {
    bad <- c(bad, "parsed/NULL mismatch")
  } else if (!is.null(expect_hours) && !isTRUE(all.equal(as.numeric(h), expect_hours))) {
    bad <- c(bad, sprintf("hours = %s, expected %s",
                          paste(h, collapse = ","), paste(expect_hours, collapse = ",")))
  }
  if (!identical(u, expect_uneven)) bad <- c(bad, sprintf("uneven = %s", u))
  if (!is.null(expect_cum) && !isTRUE(all.equal(as.numeric(cum), expect_cum)))
    bad <- c(bad, sprintf("cumulative = %s, expected %s",
                          paste(cum, collapse = ","), paste(expect_cum, collapse = ",")))
  if (length(bad)) {
    failures <<- failures + 1
    cat(sprintf("FAIL: %-10s %s\n", label, paste(bad, collapse = "; ")))
  } else {
    cat(sprintf("ok  : %-10s %s\n", label,
                if (is.null(h)) "correctly refused" else paste(h, collapse = ", ")))
  }
}

# --- shapes the column names actually take -----------------------------------
chk("hours",    c("Base9h", "Base10h", "Base11h"),               c(9, 10, 11),           FALSE)
chk("h+min",    c("Base7h30", "Base8h30", "Base9h30"),           c(7.5, 8.5, 9.5),       FALSE)
chk("dutch u",  c("KSS_9u_dag1", "KSS_11u_dag1", "KSS_13u_dag1"), c(9, 11, 13),          FALSE)
chk("HH:MM",    sprintf("%02d:00", 0:5),                          0:5,                   FALSE)
chk("decimal",  c("X8.0", "X8.25", "X8.5", "X8.75"),              c(8, 8.25, 8.5, 8.75), FALSE)

# --- the case the whole feature exists for: unequal spacing across midnight ---
chk("uneven",   c("KSS_20u", "KSS_22u", "KSS_2u", "KSS_8u"),      c(20, 22, 2, 8), TRUE,
    expect_cum = c(0, 2, 6, 12))

# --- shapes that must be REFUSED rather than mistaken for hours ---------------
# T1..T30 parses to 1..30 under the "any number" fallback; treating that as
# clock time is exactly the bug this guard exists to prevent.
chk("T1..T30",  paste0("T", 1:30),                                NULL, FALSE)
chk("no digits", c("pre", "post", "followup"),                    NULL, FALSE)
chk("duplicate", c("Base9h", "Base9h", "Base10h"),                c(9, 9, 10), FALSE)
if (!is.null(fck_cumulative_hours(c("Base9h", "Base9h", "Base10h")))) {
  failures <- failures + 1
  cat("FAIL: duplicate timestamps must not yield a time axis (zero-width step)\n")
} else {
  cat("ok  : duplicate  zero-width step refused\n")
}

if (failures) { cat("\n", failures, " failure(s).\n", sep = ""); quit(status = 1) }
cat("\nClock helper tests passed.\n")
