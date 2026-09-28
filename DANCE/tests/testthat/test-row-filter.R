# ==============================================================================
# tests/testthat/test-row-filter.R
#
# The participant/group selection on the Data Import tab
# (server/02d_helpers_rowfilter.R). What must hold:
#   * a row is analysed unless its participant, or its level of any grouping
#     variable, was excluded;
#   * every vector that travels with the curves is cut the same way, and an
#     excluded group does not survive as an empty factor level;
#   * exclusions are keyed by participant and level, so rebuilding the frame
#     (a re-confirm that drops or reorders rows) re-applies them to the SAME
#     people -- and never to somebody else;
#   * a selection that leaves fewer than two curves is refused, not applied.
# ==============================================================================
`%||%` <- function(a, b) if (is.null(a)) b else a
app_dir <- if (dir.exists("server")) "." else if (dir.exists("../../server")) "../.." else "DANCE"
source(file.path(app_dir, "server/02d_helpers_rowfilter.R"))

# 8 curves from 6 participants (P2 and P5 contribute two sessions each), a
# 3-level group with one unlabelled curve, a sex variable, and an id-like code
# that is distinct on every row.
mk <- function(ids = TRUE) {
  set.seed(4)
  n <- 8
  dat <- matrix(rnorm(n * 6), n, 6, dimnames = list(NULL, sprintf("%02d:00", 0:5)))
  cov <- data.frame(
    ID    = c("P1", "P2", "P2", "P3", "P4", "P5", "P5", "P6"),
    Group = factor(c("A", "B", "B", "A", "C", "B", "B", NA), levels = c("A", "B", "C")),
    Sex   = c("F", "M", "M", "F", "M", "F", "F", "M"),
    Code  = paste0("x", 1:8),
    stringsAsFactors = FALSE)
  gv <- cov; for (nm in names(gv)) gv[[nm]] <- as.factor(gv[[nm]])
  dance_import_frame(dat, subject_ids = if (ids) cov$ID else NULL, covariates = cov,
                     group_variables = gv, group_labels = cov$Group,
                     row_index = c(1, 2, 3, 5, 6, 7, 8, 9),   # file row 4 was dropped at import
                     time_labels = colnames(dat), time_numeric = 1:6, time_clock = 0:5)
}

test_that("participants are keyed by identifier, else by file row", {
  fr <- mk()
  expect_equal(fr$key[1:3], c("id:P1", "id:P2", "id:P2"))
  expect_equal(fr$label[2], "P2")
  fr0 <- mk(ids = FALSE)
  expect_equal(fr0$key[4], "row:5")          # the FILE row, not the position
  expect_equal(fr0$label[4], "Row 5")
  expect_equal(rownames(fr0$data)[4], "5")   # so tabs falling back to row names agree
  # a missing identifier falls back to the row for that row only
  cov_ids <- c("P1", NA, "", "P3", "P4", "P5", "P5", "P6")
  fr2 <- dance_import_frame(fr$data, subject_ids = cov_ids, row_index = fr$row_index)
  expect_equal(fr2$key[2:3], c("row:2", "row:3"))
  expect_equal(fr2$key[1], "id:P1")
})

test_that("the group picker offers the categorical variables and nothing else", {
  fr <- mk()
  v <- dance_filter_vars(fr)
  expect_setequal(names(v), c("Group", "Sex"))
  # ID is the participant id, Code is distinct on every row
  expect_false(any(c("ID", "Code") %in% names(v)))
  # factor levels keep their order, counts are of curves, missing is a level
  expect_equal(v$Group$level, c("A", "B", "C", DANCE_FILTER_MISSING))
  expect_equal(v$Group$n, c(2L, 4L, 1L, 1L))
  expect_equal(sum(v$Sex$n), 8L)
  # a numeric variable with few values counts as categorical, sorted numerically
  fr$covariates$Visit <- c(10, 2, 2, 10, 1, 2, 1, 10)
  v2 <- dance_filter_vars(fr)
  expect_equal(v2$Visit$level, c("1", "2", "10"))
})

test_that("a continuous variable is not offered as a group", {
  set.seed(9)
  n <- 40
  fr <- dance_import_frame(matrix(rnorm(n * 4), n, 4),
                           covariates = data.frame(Age = round(rnorm(n, 45, 12)),
                                                   Arm = rep(c("drug", "placebo"), 20)))
  v <- dance_filter_vars(fr)
  expect_equal(names(v), "Arm")
  # the same rule as the cosinor design picker: at most 12 distinct values
  expect_gt(length(unique(fr$covariates$Age)), 12)
})

test_that("participants are listed once, with their curves and group", {
  fr <- mk()
  p <- dance_filter_participants(fr)
  expect_equal(p$key, paste0("id:P", 1:6))
  expect_equal(p$n_rows, c(1L, 2L, 1L, 1L, 2L, 1L))
  expect_equal(p$group[6], DANCE_FILTER_MISSING)
  txt <- dance_filter_participant_text(p)
  expect_true(grepl("P2", txt[2]) && grepl("2 curves", txt[2]) && grepl("B", txt[2]))
})

test_that("ticks become a keyed exclusion record, and back", {
  fr <- mk()
  v <- dance_filter_vars(fr)
  all_lv <- unlist(lapply(names(v), function(nm) dance_filter_code(nm, v[[nm]]$level)))
  all_p <- dance_filter_participants(fr)$key
  none <- dance_filter_from_ticks(fr, all_lv, all_p)
  expect_true(dance_filter_is_empty(none))
  rf <- dance_filter_from_ticks(fr, setdiff(all_lv, "Group::C"), setdiff(all_p, "id:P3"))
  expect_equal(rf$excluded_levels, list(Group = "C"))
  expect_equal(rf$excluded_ids, "id:P3")
  # no level picker on screen: nothing is excluded by level, whatever the input
  rf2 <- dance_filter_from_ticks(fr, NULL, all_p, has_level_ui = FALSE)
  expect_true(dance_filter_is_empty(rf2))
})

test_that("a row goes when its participant OR any of its levels is excluded", {
  fr <- mk()
  fk <- dance_filter_keep(fr, excluded_ids = "id:P2",
                          excluded_levels = list(Group = c("C", DANCE_FILTER_MISSING)))
  # P2's two sessions, the C curve, the unlabelled curve
  expect_equal(which(!fk$keep), c(2L, 3L, 5L, 8L))
  expect_equal(fk$n_kept, 4L)
  expect_equal(sum(fk$by_id), 2L)
  expect_equal(fk$level_hits$Group, 2L)
  # two variables at once: a curve must pass both
  fk2 <- dance_filter_keep(fr, excluded_levels = list(Group = "A", Sex = "M"))
  expect_equal(which(fk2$keep), c(6L, 7L))   # the two B/F sessions of P5
})

test_that("applying cuts every parallel vector and drops emptied levels", {
  fr <- mk()
  fk <- dance_filter_keep(fr, excluded_levels = list(Group = c("C", DANCE_FILTER_MISSING)))
  res <- dance_filter_apply(fr, fk$keep)
  expect_true(res$ok)
  expect_equal(nrow(res$data), 6L)
  expect_equal(res$row_index, c(1L, 2L, 3L, 5L, 7L, 8L))
  expect_equal(res$subject_ids, c("P1", "P2", "P2", "P3", "P5", "P5"))
  expect_equal(nrow(res$covariates), 6L)
  expect_equal(nrow(res$group_variables), 6L)
  # the excluded group is gone from the factor, not left behind empty
  expect_equal(levels(res$group_labels), c("A", "B"))
  expect_equal(levels(res$covariates$Group), c("A", "B"))
  expect_equal(levels(res$group_variables$Group), c("A", "B"))
  # the curves are the right ones, not the first six
  expect_equal(unname(res$data), unname(fr$data[fk$keep, ]))
  expect_equal(rownames(res$data), as.character(res$row_index))
  expect_false(res$cols_dropped)
  expect_equal(res$time_clock, 0:5)
})

test_that("a time point nobody left was measured at is dropped with its label", {
  fr <- mk()
  fr$data[c(5, 8), 3] <- 5                 # only the C and the unlabelled curve
  fr$data[-c(5, 8), 3] <- NA
  fk <- dance_filter_keep(fr, excluded_levels = list(Group = c("C", DANCE_FILTER_MISSING)))
  res <- dance_filter_apply(fr, fk$keep)
  expect_true(res$cols_dropped)
  expect_equal(ncol(res$data), 5L)
  expect_equal(res$time_labels, c("00:00", "01:00", "03:00", "04:00", "05:00"))
  expect_equal(res$time_numeric, c(1L, 2L, 4L, 5L, 6L))
  expect_null(res$time_clock)              # the caller re-derives it from the labels
})

test_that("keeping everyone is the identity", {
  fr <- mk()
  res <- dance_filter_apply(fr, dance_filter_keep(fr)$keep)
  expect_true(res$ok)
  expect_equal(res$data, fr$data)
  expect_equal(res$covariates$ID, fr$covariates$ID)
  expect_equal(res$row_index, fr$row_index)
})

test_that("a selection leaving fewer than two curves is refused", {
  fr <- mk()
  fk <- dance_filter_keep(fr, excluded_levels = list(Group = c("A", "B", DANCE_FILTER_MISSING)))
  expect_equal(fk$n_kept, 1L)
  res <- dance_filter_apply(fr, fk$keep)
  expect_false(res$ok)
  expect_match(res$message, "at least two")
})

test_that("exclusions follow the PEOPLE when the frame is rebuilt", {
  fr <- mk()
  rf <- list(excluded_ids = c("id:P3", "id:P9"), excluded_levels = list(Group = "C", Diet = "vegan"))
  # a re-confirm that dropped P1's row (too few measurements) and so shifted
  # every position by one: positions would now point at the wrong people
  keep_rows <- -1
  fr2 <- dance_import_frame(fr$data[keep_rows, ], fr$subject_ids[keep_rows],
                            fr$covariates[keep_rows, ], fr$group_variables[keep_rows, ],
                            fr$group_labels[keep_rows], fr$row_index[keep_rows])
  fk2 <- dance_filter_keep(fr2, rf$excluded_ids, rf$excluded_levels)
  gone <- fr2$subject_ids[!fk2$keep]
  expect_setequal(gone, c("P3", "P4"))      # P3 by id, P4 is the only C
  # what no longer matches is reported, not silently forgotten ...
  expect_equal(fk2$stale_ids, "id:P9")
  expect_equal(fk2$stale_levels, list(Diet = "vegan"))
  expect_equal(dance_filter_n_stale(fk2), 2L)
  # ... and pruned from the stored record
  pr <- dance_filter_prune(rf, fk2)
  expect_equal(pr$excluded_ids, "id:P3")
  expect_equal(pr$excluded_levels, list(Group = "C"))
})

test_that("the description counts curves and participants", {
  fr <- mk()
  rf <- list(excluded_ids = "id:P2", excluded_levels = list(Group = "C"))
  fk <- dance_filter_keep(fr, rf$excluded_ids, rf$excluded_levels)
  d <- dance_filter_describe(fr, fk, rf)
  expect_equal(d[1], "5 of 8 curves analysed from 4 of 6 participants.")
  expect_true(any(grepl("Group = C: 1 curve", d, fixed = TRUE)))
  expect_true(any(grepl("excluded individually: 2 curves (P2)", d, fixed = TRUE)))
  # everyone kept: one line
  expect_length(dance_filter_describe(fr, dance_filter_keep(fr), NULL), 1L)
})

test_that("a per-curve vector goes back onto the file's rows", {
  x <- factor(c("k1", "k2", "k1"))
  out <- dance_to_file_rows(x, row_index = c(2L, 3L, 6L), n_file = 7L)
  expect_equal(length(out), 7L)
  expect_equal(as.character(out[c(2, 3, 6)]), c("k1", "k2", "k1"))
  expect_true(all(is.na(out[c(1, 4, 5, 7)])))
  expect_equal(levels(out), c("k1", "k2"))
  # without a record, only an exact length is trusted
  expect_equal(dance_to_file_rows(1:7, NULL, 7L), 1:7)
  expect_null(dance_to_file_rows(1:3, NULL, 7L))
  expect_null(dance_to_file_rows(1:3, c(1L, 2L, 9L), 7L))
})
