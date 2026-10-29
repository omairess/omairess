# ==============================================================================
# tests/testthat/test-excel-time-names.R
#
# Header cells formatted as times (8:00, 10:00, ..., 0:00, ..., 16:00 the next
# day) reach the app through readxl as fractions of a day, with "...N" on the
# times that recur. They are read back as the clock times they display, and the
# hour parser reads HH:MM with its minutes.
# ==============================================================================
.libPaths(c("~/Rlib", .libPaths()))
app_dir <- if (dir.exists("server")) "." else if (dir.exists("../../server")) "../.." else "DANCE"
`%||%` <- function(a, b) if (is.null(a)) b else a
e <- new.env()
for (f in c("01_helpers_time", "03_helpers_clock"))
  sys.source(file.path(app_dir, "server", paste0(f, ".R")), envir = e)

# exactly what readxl returned for the user's header row
xl <- c("subject", "0.33333333333333331...5", "0.41666666666666669...6", "0.5...7",
        "0.58333333333333337", "0.66666666666666663...9", "0.75", "0.83333333333333337",
        "0.91666666666666663", "0", "8.3333333333333329E-2", "0.16666666666666666", "0.25",
        "0.33333333333333331...17", "0.41666666666666669...18", "0.5...19",
        "0.66666666666666663...20")

test_that("Excel day fractions become the clock times they display", {
  r <- e$dance_dayfraction_names(xl)
  expect_identical(r$names[1], "subject")
  expect_false(r$converted[1])
  expect_identical(r$names[-1], c("08:00", "10:00", "12:00", "14:00", "16:00", "18:00", "20:00",
                                  "22:00", "00:00", "02:00", "04:00", "06:00", "08:00 #2",
                                  "10:00 #2", "12:00 #2", "16:00 #2"))
  expect_false(anyDuplicated(r$names) > 0)
  # and they are clock times the app accepts, crossing midnight once
  h <- e$dance_clock_hours(r$names[-1])
  expect_equal(h, c(8, 10, 12, 14, 16, 18, 20, 22, 0, 2, 4, 6, 8, 10, 12, 16))
  expect_equal(e$dance_cumulative_hours(r$names[-1]),
               c(0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 32))
})

test_that("typed decimal hours and ordinary names are left alone", {
  for (nm in list(c("0.25", "0.5", "0.75"), c("1", "2", "3"), c("8.25", "9.5"),
                  c("subject", "age"), c("0.33333333333333331")))
    expect_identical(e$dance_dayfraction_names(nm)$names, nm, info = paste(nm, collapse = ","))
})

test_that("HH:MM is read with its minutes", {
  expect_equal(e$extract_hour_from_colname("08:00"), 8)
  expect_equal(e$extract_hour_from_colname("8:30"), 8.5)
  expect_equal(e$extract_hour_from_colname("MT_21:30_SD"), 21.5)
  expect_equal(e$extract_hour_from_colname("16:00 #2"), 16)
  # what it read before is unchanged
  expect_equal(e$extract_hour_from_colname("Base7h30"), 7.5)
  expect_equal(e$extract_hour_from_colname("KSS_9u_dag1"), 9)
  expect_equal(e$extract_hour_from_colname("X8.25"), 8.25)
})
