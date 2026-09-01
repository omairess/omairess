# ==============================================================================
# tests/missing_status_test.R — the observed / interpolated / extrapolated split
#
# This is the classification the whole missing-data panel rests on, and getting
# it wrong is invisible: mislabel an extrapolated point as interpolated and the
# map simply reassures you about a value nothing constrains. The cases below are
# the ones that decide it, including a leading/trailing-gap pattern taken from a
# real staggered sleep-deprivation protocol.
#
# Run with:   Rscript tests/missing_status_test.R      (from the FCK directory)
# Needs no packages.
# ==============================================================================

app_dir <- if (dir.exists("server")) "." else "FCK"
eval(parse(file.path(app_dir, "server/01_helpers_time.R"), encoding = "UTF-8"))
eval(parse(file.path(app_dir, "server/03_helpers_clock.R"), encoding = "UTF-8"))
eval(parse(file.path(app_dir, "server/05_helpers_missing.R"), encoding = "UTF-8"))

failures <- 0
O <- FCK_FILL_OBSERVED; I <- FCK_FILL_INTERPOLATED; E <- FCK_FILL_EXTRAPOLATED

chk <- function(label, row, expect) {
  got <- as.vector(fck_fill_status(matrix(row, nrow = 1)))
  if (!identical(as.integer(got), as.integer(expect))) {
    failures <<- failures + 1
    cat(sprintf("FAIL: %-22s got %s, expected %s\n", label,
                paste(got, collapse = ""), paste(expect, collapse = "")))
  } else {
    cat(sprintf("ok  : %-22s %s\n", label,
                paste(c("O", "I", "E")[got + 1L], collapse = "")))
  }
}

chk("complete row",        c(1, 2, 3, 4),           c(O, O, O, O))
chk("interior gap",        c(1, NA, NA, 4),         c(O, I, I, O))
chk("leading gap",         c(NA, NA, 3, 4),         c(E, E, O, O))
chk("trailing gap",        c(1, 2, NA, NA),         c(O, O, E, E))
chk("both ends",           c(NA, 2, NA, 4, NA),     c(E, O, I, O, E))
chk("single observation",  c(NA, 2, NA),            c(E, O, E))
chk("nothing observed",    c(NA, NA, NA),           c(E, E, E))

# the shape that motivated the feature: a subject who joins late and leaves
# early in a long staggered protocol — most of the row is beyond their data
chk("staggered protocol",
    c(NA, NA, 3, 4, NA, 6, NA, NA, NA, NA),
    c(E, E, O, O, I, O, E, E, E, E))

# --- per-subject summary -----------------------------------------------------
m <- rbind(c(NA, NA, 3, 4, NA, 6, NA, NA, NA, NA),
           c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10))
st  <- fck_fill_status(m)
tab <- fck_fill_per_subject(st, hours_per_step = 2)

expect <- function(what, got, want) {
  if (!isTRUE(all.equal(got, want))) {
    failures <<- failures + 1
    cat(sprintf("FAIL: %-22s got %s, expected %s\n", what, format(got), format(want)))
  } else cat(sprintf("ok  : %-22s %s\n", what, format(got)))
}
expect("observed count",     tab$Observed[1],              3L)
expect("interpolated count", tab$Interpolated[1],          1L)
expect("extrapolated count", tab$Extrapolated[1],          6L)
expect("hours extrapolated", tab$`Hours extrapolated`[1],  12)
expect("longest interior gap", tab$`Longest gap`[1],       1L)
expect("first observed idx", tab$`First observed`[1],      3L)
expect("last observed idx",  tab$`Last observed`[1],       6L)
expect("complete row filled %", tab$`% filled`[2],         0)

# a gap must never be counted as extrapolation, nor the reverse
if (sum(st == FCK_FILL_INTERPOLATED) + sum(st == FCK_FILL_EXTRAPOLATED) !=
    sum(is.na(m))) {
  failures <- failures + 1
  cat("FAIL: filled cells do not add up to the missing cells\n")
} else cat("ok  : filled cells add up to the missing cells\n")

if (failures) { cat("\n", failures, " failure(s).\n", sep = ""); quit(status = 1) }
cat("\nMissing-data status tests passed.\n")
