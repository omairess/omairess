# ==============================================================================
# server/02d_helpers_rowfilter.R — which participants and groups are analysed
# ==============================================================================
# The Data Import tab builds ONE analysis frame: the curve matrix and every
# vector that travels with it row for row -- participant ids, the scalar
# variables, their factor copies, the primary grouping factor, and the row of
# the imported file each curve came from. Every analysis reads that frame
# through values$data and friends.
#
# This file decides which ROWS of the frame the analyses see. It is pure -- a
# frame in, a frame out -- so the import tab, a restored session and the tests
# apply one rule:
#
#   a row is analysed unless the analyst excluded its PARTICIPANT, or excluded
#   its LEVEL of any grouping variable.
#
# Exclusions are stored rather than inclusions, and they are keyed by
# participant and by (variable, level), never by row position. Re-confirming
# the variable selection rebuilds the frame -- and can drop rows -- so a
# positional record would silently exclude somebody else; a keyed record
# re-applies to the same people, and a participant or level that is new to the
# frame is included by default. Nothing is recomputed: excluded rows are simply
# not handed to the analyses, and every analysis is re-run on what is left.
# ==============================================================================

# The label a row with no value for a grouping variable is listed under, so an
# unlabelled row can be excluded like any other level.
DANCE_FILTER_MISSING <- "(missing)"

# Separator between a variable and a level in the group picker's values. A
# column name containing it would not decode, which is why decoding splits on
# the FIRST occurrence: a level may contain it, a variable name may not.
DANCE_FILTER_SEP <- "::"

# ---- the frame --------------------------------------------------------------
# Snapshotted once per import (and per sample dataset), so a selection can
# always be undone without reading the file again.
#
# `row_index` is the row of the imported (subjects-in-rows) frame each curve
# came from, carried through every row the import step drops. It is what keys a
# participant when the file has no identifier, and what lets a module that
# reads the raw frame (the repeated-measures pickers, the cluster-membership
# export) find the right row after rows were dropped or excluded.
#
# `file_rows = FALSE` says the positions are NOT known to be rows of the
# imported frame (a session saved before this record existed, after the import
# step had dropped rows). The frame still works -- a curve is then named by its
# position -- but it hands no row record on, so nothing reading the raw frame
# can be pointed at the wrong row by it.
dance_import_frame <- function(data, subject_ids = NULL, covariates = NULL,
                               group_variables = NULL, group_labels = NULL,
                               row_index = NULL, time_labels = NULL,
                               time_numeric = NULL, time_clock = NULL,
                               file_rows = TRUE) {
  data <- as.matrix(data)
  n <- nrow(data)
  if (is.null(row_index) || length(row_index) != n) row_index <- seq_len(n)
  row_index <- as.integer(row_index)
  # a curve is always nameable by its file row, so the tabs that fall back to
  # row names cannot renumber it after an exclusion
  if (is.null(rownames(data))) rownames(data) <- as.character(row_index)
  ids <- if (!is.null(subject_ids) && length(subject_ids) == n) as.character(subject_ids) else NULL
  has_id <- if (is.null(ids)) rep(FALSE, n) else !(is.na(ids) | !nzchar(trimws(ids)))
  key <- paste0("row:", row_index)
  key[has_id] <- paste0("id:", ids[has_id])
  # what a person reads in the picker: the identifier, else the name the file
  # gave the row (a transposed "long" file names its subjects), else its number
  rn <- rownames(data)
  own_name <- !identical(rn, as.character(row_index))
  label <- if (own_name) ifelse(grepl("^[0-9]+$", rn), paste("Row", rn), rn)
           else paste(if (isTRUE(file_rows)) "Row" else "Curve", row_index)
  label[has_id] <- ids[has_id]
  list(data = data, subject_ids = ids, covariates = covariates,
       group_variables = group_variables, group_labels = group_labels,
       row_index = row_index, key = key, label = label,
       time_labels = time_labels, time_numeric = time_numeric,
       time_clock = time_clock, file_rows = isTRUE(file_rows))
}

# One column of the frame's scalar variables, original type first.
dance_filter_column <- function(frame, nm) {
  for (src in list(frame$covariates, frame$group_variables))
    if (!is.null(src) && nm %in% names(src)) return(src[[nm]])
  NULL
}

dance_filter_labels <- function(x) {
  lab <- as.character(x)
  lab[is.na(lab) | !nzchar(trimws(lab))] <- DANCE_FILTER_MISSING
  lab
}

# ---- what a selection can be made on ---------------------------------------
# The categorical scalar variables: factor, character or logical columns, or
# numeric ones with at most 12 distinct values -- the rule the cosinor design
# picker already uses, so a variable is a "group" in one place exactly when it
# is in the other. A variable is offered only when it actually splits the rows
# (two or more levels, counting a missing label as a level of its own), when it
# has at most `max_levels` levels, and when it is not an identifier: a column
# that is the participant id, or that is distinct on every row, is not a
# grouping of anything.
dance_filter_vars <- function(frame, max_levels = 30L) {
  nms <- unique(c(names(frame$covariates), names(frame$group_variables)))
  n <- nrow(frame$data)
  out <- list()
  for (nm in nms) {
    x <- dance_filter_column(frame, nm)
    if (is.null(x) || length(x) != n) next
    lab <- dance_filter_labels(x)
    u <- unique(lab)
    obs <- x[!is.na(x)]
    cat_like <- is.factor(x) || is.character(x) || is.logical(x) ||
      length(unique(obs)) <= 12L
    if (!cat_like || length(u) < 2L || length(u) > max_levels) next
    if (!is.null(frame$subject_ids) && identical(as.character(x), frame$subject_ids)) next
    if (n > 2L && length(u) == n) next
    real <- setdiff(u, DANCE_FILTER_MISSING)
    lv <- if (is.factor(x)) levels(x)[levels(x) %in% real]
          else if (is.numeric(x)) as.character(sort(unique(obs)))
          else sort(real)
    lv <- c(lv[lv %in% real], if (DANCE_FILTER_MISSING %in% u) DANCE_FILTER_MISSING)
    out[[nm]] <- data.frame(level = lv,
                            n = vapply(lv, function(l) sum(lab == l), integer(1)),
                            stringsAsFactors = FALSE, row.names = NULL)
  }
  out
}

# ---- who can be selected ----------------------------------------------------
# One entry per participant (or per row when the file has no identifier), in
# file order, with how many curves they contribute and the primary group they
# carry -- all of it shown in the picker, because "S017" alone does not tell
# the analyst whether that is the person they meant.
dance_filter_participants <- function(frame) {
  k <- frame$key
  first <- !duplicated(k)
  n_rows <- as.integer(table(factor(k, levels = unique(k))))
  grp <- if (!is.null(frame$group_labels) && length(frame$group_labels) == length(k))
    vapply(split(dance_filter_labels(frame$group_labels), factor(k, levels = unique(k))),
           function(g) paste(unique(g), collapse = "/"), character(1))
  else rep(NA_character_, sum(first))
  data.frame(key = k[first], label = frame$label[first], n_rows = n_rows,
             group = unname(grp), stringsAsFactors = FALSE, row.names = NULL)
}

# The picker text for each participant.
dance_filter_participant_text <- function(parts) {
  paste0(parts$label,
         ifelse(is.na(parts$group), "", paste0(" | ", parts$group)),
         ifelse(parts$n_rows > 1L, paste0(" (", parts$n_rows, " curves)"), ""))
}

# ---- ticks <-> exclusions ---------------------------------------------------
dance_filter_code <- function(var, level) {
  # paste0() turns a zero-length argument into "", which would make a code for
  # no level at all
  if (!length(level)) return(character(0))
  paste0(var, DANCE_FILTER_SEP, level)
}

# The exclusion record implied by what is ticked in the two pickers. A picker
# that is not on screen (no categorical variable to offer) excludes nothing;
# one that is on screen with nothing ticked excludes everything it lists.
dance_filter_from_ticks <- function(frame, ticked_levels = NULL, ticked_keys = NULL,
                                    has_level_ui = TRUE, has_part_ui = TRUE) {
  excluded_levels <- list()
  if (isTRUE(has_level_ui)) {
    vars <- dance_filter_vars(frame)
    for (nm in names(vars)) {
      lv <- vars[[nm]]$level
      ex <- lv[!dance_filter_code(nm, lv) %in% ticked_levels]
      if (length(ex)) excluded_levels[[nm]] <- ex
    }
  }
  excluded_ids <- if (isTRUE(has_part_ui))
    setdiff(dance_filter_participants(frame)$key, ticked_keys) else character(0)
  list(excluded_ids = excluded_ids, excluded_levels = excluded_levels)
}

dance_filter_is_empty <- function(rf) {
  is.null(rf) || (!length(rf$excluded_ids) &&
                    !length(unlist(rf$excluded_levels, use.names = FALSE)))
}

# ---- the rule ---------------------------------------------------------------
# Which rows are kept, and why each dropped row was dropped. Exclusions that
# name a participant or a level no longer in the frame are returned as `stale`
# rather than silently forgotten, so the caller can say so.
dance_filter_keep <- function(frame, excluded_ids = character(0),
                              excluded_levels = list()) {
  n <- nrow(frame$data)
  by_id <- frame$key %in% excluded_ids
  by_level <- rep(FALSE, n)
  level_hits <- list()
  stale_levels <- list()
  for (nm in names(excluded_levels)) {
    ex <- excluded_levels[[nm]]
    if (!length(ex)) next
    x <- dance_filter_column(frame, nm)
    if (is.null(x) || length(x) != n) { stale_levels[[nm]] <- ex; next }
    lab <- dance_filter_labels(x)
    hit <- lab %in% ex
    by_level <- by_level | hit
    level_hits[[nm]] <- sum(hit)
    gone <- setdiff(ex, unique(lab))
    if (length(gone)) stale_levels[[nm]] <- gone
  }
  keep <- !(by_id | by_level)
  list(keep = keep, by_id = by_id, by_level = by_level, level_hits = level_hits,
       stale_ids = setdiff(excluded_ids, frame$key), stale_levels = stale_levels,
       n_total = n, n_kept = sum(keep))
}

# ---- applying it ------------------------------------------------------------
# The kept rows of every parallel vector. Factor columns lose the levels nobody
# kept carries, so an excluded group does not survive as an empty level -- an
# empty level is a group with no curves, and the group tests downstream would
# either fail on it or count it. A time point that none of the remaining
# curves was measured at is dropped, exactly as the import step drops one that
# nobody was measured at; the caller re-derives the clock times from the
# remaining labels.
dance_filter_apply <- function(frame, keep) {
  keep <- as.logical(keep)
  n_kept <- sum(keep)
  # Refusals are about what a SELECTION does. Keeping everyone is the identity,
  # whatever the frame holds; judging the frame itself is the import step's job.
  excluding <- !all(keep)
  if (excluding && n_kept < 2L)
    return(list(ok = FALSE, n_kept = n_kept, message = sprintf(
      "That selection leaves %d curve%s; at least two are needed for any analysis.",
      n_kept, if (n_kept == 1L) "" else "s")))
  sub_df <- function(df) {
    if (is.null(df)) return(NULL)
    df <- df[keep, , drop = FALSE]
    for (nm in names(df)) if (is.factor(df[[nm]])) df[[nm]] <- droplevels(df[[nm]])
    df
  }
  dat <- frame$data[keep, , drop = FALSE]
  cols <- if (excluding) colSums(!is.na(dat)) > 0 else rep(TRUE, ncol(dat))
  if (excluding && sum(cols) < 2L)
    return(list(ok = FALSE, n_kept = n_kept, message = paste(
      "The curves that selection keeps share fewer than two measured time points;",
      "nothing can be analysed on them.")))
  sub_vec <- function(v) if (!is.null(v) && length(v) == length(cols)) v[cols] else v
  gl <- frame$group_labels
  if (!is.null(gl) && length(gl) == length(keep)) {
    gl <- gl[keep]
    if (is.factor(gl)) gl <- droplevels(gl)
  }
  list(ok = TRUE, n_kept = n_kept,
       data = dat[, cols, drop = FALSE],
       subject_ids = if (is.null(frame$subject_ids)) NULL else frame$subject_ids[keep],
       covariates = sub_df(frame$covariates),
       group_variables = sub_df(frame$group_variables),
       group_labels = gl,
       row_index = frame$row_index[keep],
       cols_keep = cols, cols_dropped = !all(cols),
       time_labels = sub_vec(frame$time_labels),
       time_numeric = sub_vec(frame$time_numeric),
       time_clock = if (all(cols)) frame$time_clock else NULL)
}

# ---- saying what was done ---------------------------------------------------
# Plain lines for the status panel, the session note and the report. Counts are
# of CURVES (rows); participants are counted separately when the file has an
# identifier, because one person can contribute several curves.
dance_filter_describe <- function(frame, fk, rf = NULL) {
  ids <- frame$subject_ids
  who <- if (!is.null(ids))
    sprintf(" from %d of %d participants",
            length(unique(ids[fk$keep & !is.na(ids)])), length(unique(ids[!is.na(ids)])))
  else ""
  out <- sprintf("%d of %d curves analysed%s.", fk$n_kept, fk$n_total, who)
  if (fk$n_kept == fk$n_total) return(out)
  lv <- if (is.null(rf$excluded_levels)) list() else rf$excluded_levels
  for (nm in names(fk$level_hits))
    out <- c(out, sprintf("  excluded by %s = %s: %d curve%s", nm,
                          paste(lv[[nm]], collapse = ", "), fk$level_hits[[nm]],
                          if (fk$level_hits[[nm]] == 1L) "" else "s"))
  n_id <- sum(fk$by_id)
  if (n_id) {
    lab <- unique(frame$label[fk$by_id])
    out <- c(out, sprintf("  excluded individually: %d curve%s (%s%s)", n_id,
                          if (n_id == 1L) "" else "s",
                          paste(utils::head(lab, 12), collapse = ", "),
                          if (length(lab) > 12) sprintf(", and %d more", length(lab) - 12) else ""))
  }
  out
}

# Exclusions in a record that match nothing in the frame any more -- a
# participant or level that a re-confirmed selection no longer contains.
dance_filter_n_stale <- function(fk) length(fk$stale_ids) + length(unlist(fk$stale_levels))

# The record with its stale entries removed, so what is stored is what applies.
dance_filter_prune <- function(rf, fk) {
  if (is.null(rf)) return(NULL)
  ids <- setdiff(rf$excluded_ids, fk$stale_ids)
  lv <- list()
  for (nm in names(rf$excluded_levels)) {
    keep <- setdiff(rf$excluded_levels[[nm]], fk$stale_levels[[nm]])
    if (length(keep)) lv[[nm]] <- keep
  }
  out <- list(excluded_ids = ids, excluded_levels = lv)
  if (dance_filter_is_empty(out)) NULL else out
}

# ---- back to the file's rows ------------------------------------------------
# A per-curve vector placed on the rows of the imported frame, with NA for rows
# that are not analysed (dropped at import or excluded here). NULL when the
# alignment cannot be vouched for.
dance_to_file_rows <- function(x, row_index, n_file) {
  if (is.null(x) || is.null(n_file)) return(NULL)
  if (!is.null(row_index) && length(row_index) == length(x) &&
      all(row_index >= 1L & row_index <= n_file)) {
    out <- rep(NA, n_file)
    if (is.factor(x)) {
      out <- factor(out, levels = levels(x)); out[row_index] <- x
    } else out[row_index] <- x
    return(out)
  }
  if (length(x) == n_file) return(x)
  NULL
}
