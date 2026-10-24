# ============================================================================
# 02_diary_helpers.R — import, sleep-diary derivation and composition building
#
# Plain R, no Shiny. The app calls these functions, and the exported script
# embeds this file verbatim and calls them with the same arguments, so the
# script reproduces the app's numbers exactly (house rule 6).
#
# Sleep-diary definitions follow the Consensus Sleep Diary
#   Carney, C. E., et al. (2012). The Consensus Sleep Diary: standardizing
#   prospective sleep self-monitoring. Sleep, 35(2), 287–302.
#   doi:10.5665/sleep.1642
# Zero replacement follows the multiplicative replacement of
#   Martín-Fernández, J. A., Barceló-Vidal, C., & Pawlowsky-Glahn, V. (2003).
#   Dealing with zeros and missing values in compositional data sets using
#   nonparametric imputation. Mathematical Geology, 35(3), 253–278.
#   doi:10.1023/A:1023866030544
# ============================================================================

`%||%` <- function(a, b) {
  if (is.null(a) || length(a) == 0 || identical(a, "")) b else a
}

# ---- names -----------------------------------------------------------------

# Syntactic, accent-free column names that brms formulas accept
# ("Qualité du sommeil" -> "Qualite_du_sommeil"). Case is kept so FJ/QR/QS
# stay recognisable.
clean_names <- function(x) {
  x <- stringi::stri_trans_general(as.character(x), "Latin-ASCII")
  x <- trimws(x)
  x <- gsub("[^A-Za-z0-9]+", "_", x)
  x <- gsub("^_+|_+$", "", x)
  x[is.na(x) | x == ""] <- "V"
  x <- ifelse(grepl("^[A-Za-z]", x), x, paste0("V", x))
  make.unique(x, sep = "_")
}

# ---- reading files into character grids ------------------------------------

# Every source is first turned into a character matrix ("grid") exactly as it
# appears in the file, so one set of layout rules serves Excel, CSV and JSON.

sniff_delim <- function(path) {
  lines <- readLines(path, n = 20, warn = FALSE, encoding = "UTF-8")
  cand  <- c(";" = ";", "," = ",", "\t" = "\t", "|" = "|")
  hits  <- vapply(cand, function(d)
    sum(lengths(regmatches(lines, gregexpr(d, lines, fixed = TRUE)))),
    numeric(1))
  names(cand)[which.max(hits)]
}

df_to_grid <- function(df) {
  df <- as.data.frame(df, stringsAsFactors = FALSE, check.names = FALSE)
  vals <- lapply(df, function(col) {
    if (is.list(col)) {
      vapply(col, function(v) paste(unlist(v), collapse = "; "), character(1))
    } else if (inherits(col, "haven_labelled")) {
      as.character(haven::as_factor(col))
    } else if (inherits(col, "POSIXt")) {
      format(col, "%Y-%m-%d %H:%M:%S")
    } else {
      as.character(col)
    }
  })
  grid <- rbind(names(df), do.call(cbind, vals))
  if (is.null(dim(grid))) grid <- matrix(grid, nrow = 1)
  grid[grid %in% c("", "NA", "NaN")] <- NA
  unname(grid)
}

# All record tables inside a parsed JSON document, keyed by their path
# ("entries", "data$nights", ...). Handles arrays of objects and objects
# keyed by date ({"2025-03-01": {...}, ...}).
json_tables <- function(x, path = "") {
  out <- list()
  key <- if (path == "") "<root>" else path
  if (is.data.frame(x)) {
    out[[key]] <- jsonlite::flatten(x)
    for (nm in names(x)) {
      if (is.list(x[[nm]]) && !is.data.frame(x[[nm]])) {
        inner <- Filter(is.data.frame, x[[nm]])
        if (length(inner)) {
          sub <- data.table::rbindlist(inner, fill = TRUE, idcol = ".parent")
          out[[paste0(key, "$", nm)]] <- as.data.frame(sub)
        }
      }
    }
    return(out)
  }
  if (is.list(x) && length(x)) {
    recordish <- vapply(x, function(el)
      is.list(el) && !is.data.frame(el) && length(el) > 0 &&
        all(vapply(el, function(v) is.null(v) || length(v) <= 1, logical(1))),
      logical(1))
    if (length(x) >= 2 && all(recordish)) {
      rows <- lapply(x, function(el) {
        el[vapply(el, is.null, logical(1))] <- NA
        as.data.frame(el, stringsAsFactors = FALSE, check.names = FALSE)
      })
      tab <- data.table::rbindlist(rows, fill = TRUE)
      if (!is.null(names(x))) tab <- cbind(.key = names(x), tab)
      out[[key]] <- as.data.frame(tab)
    }
    nms <- names(x) %||% as.character(seq_along(x))
    for (i in seq_along(x)) {
      if (is.list(x[[i]])) {
        child <- if (path == "") nms[i] else paste0(path, "$", nms[i])
        out <- c(out, json_tables(x[[i]], child))
      }
    }
  }
  out
}

# Top-level scalar fields of a JSON document (participant code, name...).
json_meta <- function(x) {
  if (is.data.frame(x) || !is.list(x)) return(list())
  keep <- vapply(x, function(v) is.atomic(v) && length(v) == 1, logical(1))
  x[keep]
}

# Read one file into a list of grids (one per sheet / table).
read_grids <- function(path, name = basename(path), sheets = "first",
                       delim = "auto", json_table = NULL) {
  ext <- tolower(tools::file_ext(name))
  if (ext %in% c("xlsx", "xls", "xlsm")) {
    all_sheets <- readxl::excel_sheets(path)
    use <- if (identical(sheets, "first")) all_sheets[1]
           else if (identical(sheets, "all")) all_sheets
           else all_sheets[grepl(sheets, all_sheets, ignore.case = TRUE)]
    if (!length(use)) {
      warning(sprintf("%s: no sheet matches '%s'", name, sheets))
      return(list())
    }
    out <- lapply(use, function(s) {
      g <- suppressMessages(readxl::read_excel(
        path, sheet = s, col_names = FALSE, col_types = "text",
        .name_repair = "minimal"))
      g <- as.matrix(g)
      if (!length(g)) g <- matrix(character(0), 0, 0)
      unname(g)
    })
    names(out) <- use
    return(out)
  }
  if (ext %in% c("csv", "txt", "tsv")) {
    d <- if (identical(delim, "auto")) sniff_delim(path) else delim
    g <- readr::read_delim(path, delim = d, col_names = FALSE,
                           col_types = readr::cols(.default = "c"),
                           trim_ws = TRUE, show_col_types = FALSE,
                           name_repair = "minimal")
    out <- list(unname(as.matrix(g)))
    names(out) <- "table"
    return(out)
  }
  if (ext == "json") {
    js   <- jsonlite::fromJSON(path, simplifyVector = TRUE, flatten = TRUE)
    tabs <- json_tables(js)
    if (!length(tabs)) stop(name, ": no table of records found in the JSON.")
    pick <- if (!is.null(json_table) && json_table %in% names(tabs)) json_table
            else names(tabs)[which.max(vapply(tabs, nrow, numeric(1)))]
    tab  <- tabs[[pick]]
    meta <- json_meta(js)
    for (m in names(meta)) {
      if (!m %in% names(tab)) tab[[paste0("meta_", m)]] <- as.character(meta[[m]])
    }
    out <- list(df_to_grid(tab))
    names(out) <- pick
    return(out)
  }
  if (ext == "sav") {
    out <- list(df_to_grid(haven::read_sav(path)))
  } else if (ext == "rds") {
    out <- list(df_to_grid(readRDS(path)))
  } else if (ext %in% c("rdata", "rda")) {
    e <- new.env()
    load(path, envir = e)
    dfs <- Filter(function(o) is.data.frame(get(o, envir = e)), ls(e))
    if (!length(dfs)) stop(name, ": no data frame inside the .RData file.")
    obj <- json_table %||% dfs[1]
    out <- list(df_to_grid(get(obj, envir = e)))
  } else {
    stop("Unsupported file type: .", ext)
  }
  names(out) <- "table"
  out
}

# ---- layouts ---------------------------------------------------------------

drop_empty <- function(g) {
  if (!length(g)) return(g)
  g[!is.na(g) & trimws(g) == ""] <- NA
  keep_r <- rowSums(!is.na(g)) > 0
  keep_c <- colSums(!is.na(g)) > 0
  g[keep_r, keep_c, drop = FALSE]
}

# Rows of titles above the real header: pick the first row that is at least
# 60% as full as the fullest of the first 30 rows.
guess_header_row <- function(g) {
  n <- min(nrow(g), 30)
  if (n == 0) return(0L)
  filled <- rowSums(!is.na(g[seq_len(n), , drop = FALSE]))
  as.integer(which(filled >= 0.6 * max(filled))[1] - 1)
}

is_numberish <- function(x) {
  x <- trimws(gsub(",", ".", x))
  !is.na(suppressWarnings(as.numeric(x)))
}

# Rows = variables, columns = days (typical paper-to-Excel CBT-I diary).
items_in_rows <- function(g, skip = 0) {
  if (skip > 0) g <- g[-seq_len(skip), , drop = FALSE]
  g <- drop_empty(g)
  if (!nrow(g) || ncol(g) < 2) return(data.frame())
  # The label column is the first of the first three columns that is mostly
  # text (item numbers or dates in column 1 are skipped).
  label_col <- 1L
  for (j in seq_len(min(3, ncol(g) - 1))) {
    v <- g[, j][!is.na(g[, j])]
    if (length(v) && mean(!is_numberish(v)) >= 0.6) { label_col <- j; break }
  }
  labels <- g[, label_col]
  labels[is.na(labels)] <- "day_label"
  vals <- g[, -seq_len(label_col), drop = FALSE]
  # Several weeks stacked under each other repeat the item labels: split into
  # blocks at each repeat of the first repeated label and stack the blocks.
  key <- tolower(stringi::stri_trans_general(trimws(labels), "Latin-ASCII"))
  dup <- unique(key[duplicated(key) & key != "day_label"])
  block <- if (length(dup)) {
    anchor <- dup[1]
    b <- cumsum(key == anchor)
    b[b == 0] <- 1
    b
  } else rep(1L, length(key))
  pieces <- lapply(split(seq_along(labels), block), function(idx) {
    m <- t(vals[idx, , drop = FALSE])
    colnames(m) <- labels[idx]
    m <- m[rowSums(!is.na(m)) > 0, , drop = FALSE]
    d <- as.data.frame(m, stringsAsFactors = FALSE, check.names = FALSE)
    names(d) <- make.unique(names(d), sep = "_")
    d
  })
  pieces <- Filter(function(d) nrow(d) > 0, pieces)
  if (!length(pieces)) return(data.frame())
  as.data.frame(data.table::rbindlist(pieces, fill = TRUE, use.names = TRUE),
                check.names = FALSE)
}

# Header row + one row per record (long), or one row per person (wide).
header_table <- function(g, skip = 0) {
  if (skip > 0) g <- g[-seq_len(skip), , drop = FALSE]
  g <- drop_empty(g)
  if (nrow(g) < 2) return(data.frame())
  hdr <- g[1, ]
  hdr[is.na(hdr)] <- paste0("V", which(is.na(hdr)))
  d <- as.data.frame(g[-1, , drop = FALSE], stringsAsFactors = FALSE)
  names(d) <- make.unique(hdr, sep = "_")
  d[rowSums(!is.na(d)) > 0, , drop = FALSE]
}

# Convert character columns that are entirely numbers (comma or point
# decimals) to numeric. Times stored by Excel become day fractions, which
# parse_clock() understands.
auto_numeric <- function(d) {
  for (nm in names(d)) {
    x <- d[[nm]]
    if (!is.character(x)) next
    y <- trimws(x)
    y[y == ""] <- NA
    ok <- !is.na(y)
    if (any(ok) && all(is_numberish(y[ok]))) {
      d[[nm]] <- as.numeric(gsub(",", ".", y))
    }
  }
  d
}

fill_down <- function(x) {
  if (!length(x)) return(x)
  for (i in seq_along(x)[-1]) if (is.na(x[i])) x[i] <- x[i - 1]
  x
}

#' Import one or more diary files into one long table (one row per night).
#'
#' @param paths      file paths
#' @param names      original file names (IDs default to these, minus extension)
#' @param layout     "long" (rows = nights), "items_in_rows" (rows = diary
#'                   items, columns = days), or "wide_days" (one row per
#'                   person, day-numbered columns)
#' @param sheets     Excel: "first", "all", or a regular expression on sheet names
#' @param skip       title rows above the header; NA = guess per sheet
#' @param delim      text files: "auto" or the delimiter
#' @param json_table JSON: path of the record table (NULL = largest)
#' @param id_from    "file", "sheet" or "column"
#' @param id_column  column holding IDs when id_from = "column"
#' @param wide_pattern regex with two groups (variable, day) for wide_days
#' @param day_position "suffix" (QS_1) or "prefix" (D1_QS) for wide_days
import_diary_files <- function(paths, names = basename(paths),
                               layout = c("long", "items_in_rows", "wide_days"),
                               sheets = "first", skip = NA, delim = "auto",
                               json_table = NULL,
                               id_from = c("file", "sheet", "column"),
                               id_column = NULL,
                               wide_pattern = NULL,
                               day_position = c("suffix", "prefix")) {
  layout       <- match.arg(layout)
  id_from      <- match.arg(id_from)
  day_position <- match.arg(day_position)
  log <- character(0)
  labels <- character(0)
  parts <- list()
  for (i in seq_along(paths)) {
    grids <- read_grids(paths[i], names[i], sheets = sheets, delim = delim,
                        json_table = json_table)
    for (s in names(grids)) {
      g <- grids[[s]]
      if (!length(g)) next
      g[!is.na(g) & trimws(g) == ""] <- NA
      sk <- if (is.na(skip)) guess_header_row(drop_empty(g)) else skip
      d <- if (layout == "items_in_rows") items_in_rows(drop_empty(g), sk)
           else header_table(drop_empty(g), sk)
      if (!nrow(d)) {
        log <- c(log, sprintf("%s [%s]: nothing readable, skipped", names[i], s))
        next
      }
      orig <- names(d)
      names(d) <- clean_names(orig)
      labels[names(d)] <- orig
      d <- cbind(.file = names[i], .sheet = s, .row = seq_len(nrow(d)), d,
                 stringsAsFactors = FALSE)
      log <- c(log, sprintf("%s [%s]: %d rows x %d variables (header rows skipped: %d)",
                            names[i], s, nrow(d), ncol(d) - 3, sk))
      parts[[length(parts) + 1]] <- d
    }
  }
  if (!length(parts)) stop("No data could be read from the selected file(s).")
  out <- as.data.frame(data.table::rbindlist(parts, fill = TRUE, use.names = TRUE),
                       stringsAsFactors = FALSE)

  # IDs
  stem <- tools::file_path_sans_ext(out$.file)
  out$ID <- switch(id_from,
    file   = stem,
    sheet  = if (length(unique(stem)) > 1) paste(stem, out$.sheet, sep = "_")
             else out$.sheet,
    column = {
      if (is.null(id_column) || !id_column %in% names(out))
        stop("Choose the column that holds participant IDs.")
      fill_down(as.character(out[[id_column]]))
    })
  out <- out[!is.na(out$ID), , drop = FALSE]

  # wide -> long: one row per participant-night
  if (layout == "wide_days") {
    pat <- wide_pattern %||% if (day_position == "suffix")
      "^(.*?)_?([0-9]{1,3})$" else "^[A-Za-z]*_?([0-9]{1,3})_(.+)$"
    day_cols <- grep(pat, setdiff(names(out), c(".file", ".sheet", ".row", "ID")),
                     value = TRUE, perl = TRUE)
    if (!length(day_cols)) stop("No column names match the day pattern: ", pat)
    names_to <- if (day_position == "suffix") c(".value", "night") else c("night", ".value")
    out <- tidyr::pivot_longer(out, cols = tidyselect::all_of(day_cols),
                               names_to = names_to, names_pattern = pat,
                               values_drop_na = FALSE)
    out <- as.data.frame(out, stringsAsFactors = FALSE)
    out$night <- as.integer(out$night)
    day_vars <- unique(base::sub(pat, if (day_position == "suffix") "\\1" else "\\2",
                           day_cols, perl = TRUE))
    out <- out[rowSums(!is.na(out[, day_vars, drop = FALSE])) > 0, , drop = FALSE]
    out <- out[order(out$ID, out$night), , drop = FALSE]
    log <- c(log, sprintf("Wide -> long: %d day-numbered columns pivoted into %d participant-nights",
                          length(day_cols), nrow(out)))
  } else {
    out$night <- stats::ave(seq_len(nrow(out)), out$ID, FUN = seq_along)
  }
  out <- auto_numeric(out)
  front <- c("ID", "night", ".file", ".sheet", ".row")
  out <- out[, c(front, setdiff(names(out), front)), drop = FALSE]
  rownames(out) <- NULL
  attr(out, "labels") <- labels
  attr(out, "log") <- c(log, sprintf("Total: %d participant(s), %d rows",
                                     length(unique(out$ID)), nrow(out)))
  out
}

# ---- clock times and durations ---------------------------------------------

#' Clock time -> minutes after midnight (0–1439).
#' Understands "23:30", "23h30", "23u30", "11:30 pm", Excel day fractions
#' (0.979), Excel date-times (45291.98), 2330, and "YYYY-mm-dd HH:MM:SS".
#' @param numeric_style how plain numbers 1–24 are read: "hh.mm" (23.30 is
#'   23:30) or "decimal_hours" (23.5 is 23:30)
parse_clock <- function(x, numeric_style = c("hh.mm", "decimal_hours")) {
  numeric_style <- match.arg(numeric_style)
  if (inherits(x, "POSIXt")) {
    lt <- as.POSIXlt(x)
    return(lt$hour * 60 + lt$min)
  }
  num_rule <- function(v) {
    out <- rep(NA_real_, length(v))
    ok  <- !is.na(v)
    frac <- ok & v >= 0 & v < 1
    out[frac] <- v[frac] * 1440
    serial <- ok & v >= 20000
    out[serial] <- (v[serial] %% 1) * 1440
    hrs <- ok & v >= 1 & v < 24.6 & !serial
    if (numeric_style == "hh.mm") {
      mm <- round((v[hrs] %% 1) * 100)
      out[hrs] <- ifelse(mm < 60, floor(v[hrs]) * 60 + mm, v[hrs] * 60)
    } else {
      out[hrs] <- v[hrs] * 60
    }
    mil <- ok & v >= 100 & v <= 2400 & v == round(v) & (v %% 100) < 60
    out[mil] <- floor(v[mil] / 100) * 60 + v[mil] %% 100
    out
  }
  if (is.numeric(x)) return(round(num_rule(x)) %% 1440)
  s <- tolower(trimws(as.character(x)))
  s[s == ""] <- NA
  res <- rep(NA_real_, length(s))
  # date-time strings: keep the time part
  dt <- !is.na(s) & grepl("^\\d{4}-\\d{2}-\\d{2}[ t]\\d{1,2}:\\d{2}", s)
  s[dt] <- base::sub("^\\d{4}-\\d{2}-\\d{2}[ t]", "", s[dt])
  # "23:30", "23h30", "23u30" (a point is left to the numeric rule, where
  # numeric_style decides between 23.30 = 23:30 and 23.5 = 23:30)
  hm <- regmatches(s, regexec("^(\\d{1,2})\\s*[:hu]\\s*(\\d{2})", s))
  for (i in which(lengths(hm) == 3)) {
    res[i] <- as.numeric(hm[[i]][2]) * 60 + as.numeric(hm[[i]][3])
  }
  honly <- is.na(res) & !is.na(s) & grepl("^\\d{1,2}\\s*[hu]$", s)
  res[honly] <- as.numeric(base::sub("\\s*[hu]$", "", s[honly])) * 60
  rest <- is.na(res) & !is.na(s)
  res[rest] <- num_rule(suppressWarnings(as.numeric(gsub(",", ".", s[rest]))))
  pm <- !is.na(res) & grepl("pm", s) & res < 720
  res[pm] <- res[pm] + 720
  am12 <- !is.na(res) & grepl("am", s) & res >= 720 & res < 780
  res[am12] <- res[am12] - 720
  round(res) %% 1440
}

#' Duration -> minutes.
#' @param unit "auto", "minutes", "hours" or "day_fraction". "auto" reads
#'   "1h30"/"1:30" as hh:mm, "20 min" as minutes, ranges "10-15" as their
#'   midpoint, and a numeric column whose values are all below 1 as Excel
#'   time fractions; other numbers as minutes.
parse_duration <- function(x, unit = c("auto", "minutes", "hours", "day_fraction")) {
  unit <- match.arg(unit)
  if (inherits(x, "POSIXt")) {
    lt <- as.POSIXlt(x)
    return(lt$hour * 60 + lt$min)
  }
  s <- tolower(trimws(stringi::stri_trans_general(as.character(x), "Latin-ASCII")))
  s[s == ""] <- NA
  s <- gsub(",", ".", s)
  res <- rep(NA_real_, length(s))
  hm <- regmatches(s, regexec("^(\\d{1,2})\\s*[:hu]\\s*(\\d{1,2})?\\s*(min)?$", s))
  is_hm <- lengths(hm) > 0
  for (i in which(is_hm)) {
    mins <- if (hm[[i]][3] == "") 0 else as.numeric(hm[[i]][3])
    res[i] <- as.numeric(hm[[i]][2]) * 60 + mins
  }
  rng <- !is_hm & !is.na(s) & grepl("^\\d+(\\.\\d+)?\\s*(-|a|tot|to)\\s*\\d+(\\.\\d+)?", s)
  for (i in which(rng)) {
    nums <- as.numeric(regmatches(s[i], gregexpr("\\d+(\\.\\d+)?", s[i]))[[1]][1:2])
    res[i] <- mean(nums)
  }
  rest <- !is_hm & !rng & !is.na(s)
  num <- suppressWarnings(as.numeric(gsub("[^0-9.]", "", s[rest])))
  num[gsub("[^0-9.]", "", s[rest]) == ""] <- NA
  had_h <- grepl("(h|uur|hour)", s[rest]) & !grepl("min", s[rest])
  mult <- switch(unit,
    minutes = 1, hours = 60, day_fraction = 1440,
    auto = {
      all_num <- suppressWarnings(as.numeric(s[!is.na(s)]))
      if (length(all_num) && all(!is.na(all_num)) && all(all_num < 1) &&
          any(all_num > 0)) 1440 else 1
    })
  res[rest] <- ifelse(had_h, num * 60, num * mult)
  res
}

# Minutes from clock a to clock b, wrapping past midnight.
clock_diff <- function(a, b) (b - a) %% 1440

parse_date_any <- function(x) {
  if (inherits(x, "Date")) return(x)
  if (inherits(x, "POSIXt")) return(as.Date(x))
  if (is.numeric(x)) {
    out <- as.Date(rep(NA_real_, length(x)), origin = "1970-01-01")
    ok <- !is.na(x) & x > 20000 & x < 80000
    out[ok] <- as.Date(floor(x[ok]), origin = "1899-12-30")
    return(out)
  }
  s <- trimws(as.character(x))
  num <- suppressWarnings(as.numeric(s))
  out <- parse_date_any(num)
  # European day-first dates; ISO dates when the year comes first
  t <- gsub("[./]", "-", base::sub("[ T].*$", "", s))
  iso <- is.na(out) & grepl("^\\d{4}-\\d{1,2}-\\d{1,2}$", t)
  out[iso] <- as.Date(t[iso], format = "%Y-%m-%d")
  dmy4 <- is.na(out) & grepl("^\\d{1,2}-\\d{1,2}-\\d{4}$", t)
  out[dmy4] <- as.Date(t[dmy4], format = "%d-%m-%Y")
  dmy2 <- is.na(out) & grepl("^\\d{1,2}-\\d{1,2}-\\d{2}$", t)
  out[dmy2] <- as.Date(t[dmy2], format = "%d-%m-%y")
  out
}

# ---- recognising diary items (NL / FR / EN) --------------------------------

DIARY_FIELDS <- c(
  bed        = "Bedtime (got into bed)",
  lights_off = "Lights off / tried to sleep",
  sol        = "Sleep-onset latency (min)",
  n_wake     = "Number of awakenings",
  waso       = "Wake after sleep onset (min)",
  final_wake = "Final awakening (clock)",
  rise       = "Out of bed (clock)",
  tst        = "Total sleep time (if recorded)",
  tib        = "Time in bed (if recorded)",
  nap        = "Naps (min)",
  date       = "Date"
)

DIARY_PATTERNS <- list(
  bed        = "coucher|couche|au_lit|heure.*lit|naar_bed|bedtijd|in_bed_gegaan|ging.*bed|bed_?time|went_to_bed|got_into_bed|into_bed",
  lights_off = "extinction|eteint|lumiere|licht.*uit|lights?_?off|try.*sleep|tried.*sleep|essay.*dorm|slapen.*probe|probe.*slapen",
  sol        = "^SOL$|latenc|endorm|inslaap|in_slaap|sleep_?onset|fall_?asleep|to_fall|minutes_?to_?sleep",
  n_wake     = "nombre.*reveil|nb.*reveil|aantal.*wak|keer.*wakker|number.*awak|n_?awak|NWAK",
  waso       = "WASO|duree.*reveil|duree.*eveil|eveil.*nuit|temps.*eveill|wakker.*duur|duur.*wakker|hoelang.*wakker|time_awake|minutes_awake|awake_during",
  final_wake = "reveil_final|heure.*reveil|derniere.*reveil|reveille|wakker_geworden|ontwaak|eindelijk_wakker|final_?awak|final_?wake|last_?wake|wake_?time|woke|wake_?up",
  rise       = "lever|leve|sortie.*lit|sorti.*lit|opgestaan|opstaan|uit_bed|out_?of_?bed|rise_?time|got_?up|get_?up|arise",
  tst        = "^TST$|temps.*total.*sommeil|duree.*sommeil|totale_slaap|slaapduur|total_sleep|sleep_?duration|slept",
  tib        = "^TIB$|temps.*au_lit|tijd_in_bed|time_in_bed",
  nap        = "sieste|dutje|middagslaap|nap",
  date       = "^date|datum|^jour$|^dag$|^night$|^nuit$|^nacht$|day_label"
)

# First column whose (clean) name matches each field's pattern. The
# bookkeeping columns made by import_diary_files() are never candidates.
guess_diary_map <- function(cols) {
  cols <- setdiff(cols, c("ID", "night", ".file", ".sheet", ".row"))
  used <- character(0)
  out <- list()
  for (f in names(DIARY_PATTERNS)) {
    hit <- grep(DIARY_PATTERNS[[f]], cols, ignore.case = TRUE, value = TRUE, perl = TRUE)
    hit <- setdiff(hit, used)
    if (f == "final_wake") hit <- hit[!grepl("nombre|nb|aantal|number|duree|duur|qualit|kwalit", hit, ignore.case = TRUE)]
    if (f == "bed") hit <- hit[!grepl("lever|sortie|uit|out", hit, ignore.case = TRUE)]
    if (length(hit)) { out[[f]] <- hit[1]; used <- c(used, hit[1]) }
  }
  out
}

# ---- deriving sleep parameters ---------------------------------------------

#' Sleep parameters per night (Consensus Sleep Diary definitions).
#'
#' TIB = out of bed - into bed; TST = (final awakening - lights off) - SOL -
#' WASO, unless TST is recorded; EMA = out of bed - final awakening;
#' PRE = lights off - into bed. A night is valid when 0 < TST <= TIB < 1440.
#'
#' @param map named list: field -> column name (see DIARY_FIELDS)
#' @param waso_blank_zero treat a blank WASO cell as 0 min
derive_sleep_diary <- function(data, map, clock_numeric = "hh.mm",
                               duration_unit = "auto", waso_blank_zero = FALSE) {
  col <- function(f) {
    nm <- map[[f]]
    if (is.null(nm) || identical(nm, "") || !nm %in% names(data)) NULL else data[[nm]]
  }
  n <- nrow(data)
  clk <- function(f) { v <- col(f); if (is.null(v)) rep(NA_real_, n) else parse_clock(v, clock_numeric) }
  dur <- function(f) { v <- col(f); if (is.null(v)) rep(NA_real_, n) else parse_duration(v, duration_unit) }

  bed  <- clk("bed")
  off  <- clk("lights_off"); off[is.na(off)] <- bed[is.na(off)]
  wake <- clk("final_wake")
  rise <- clk("rise");       rise[is.na(rise)] <- wake[is.na(rise)]
  sol  <- dur("sol")
  waso <- dur("waso")
  if (isTRUE(waso_blank_zero)) waso[is.na(waso)] <- 0
  nap  <- dur("nap"); nap[is.na(nap)] <- 0

  tib <- if (!is.null(col("tib"))) dur("tib") else clock_diff(bed, rise)
  tib[!is.na(tib) & tib == 0] <- NA
  tst <- if (!is.null(col("tst"))) dur("tst")
         else clock_diff(off, wake) - sol - waso

  out <- data
  out$bed_clock   <- bed
  out$off_clock   <- off
  out$wake_clock  <- wake
  out$rise_clock  <- rise
  out$TIB  <- tib
  out$TST  <- tst
  out$SOL  <- sol
  out$WASO <- waso
  out$PRE  <- ifelse(is.na(bed) | is.na(off), 0, clock_diff(bed, off))
  out$EMA  <- ifelse(is.na(wake) | is.na(rise), 0, clock_diff(wake, rise))
  out$NAP  <- nap
  out$SE   <- round(100 * tst / tib, 1)
  if (!is.null(col("date"))) out$date <- parse_date_any(col("date"))

  why <- rep("", n)
  why[is.na(tib)] <- "time in bed missing"
  why[why == "" & is.na(tst)] <- "sleep time missing (check SOL / WASO / clocks)"
  why[why == "" & tst <= 0] <- "sleep time <= 0"
  why[why == "" & tst > tib] <- "sleep time > time in bed"
  why[why == "" & tib >= 1440] <- "time in bed >= 24 h"
  why[why == "" & out$EMA > 720] <- "out of bed before final awakening"
  why[why == "" & out$PRE > 720] <- "lights off before going to bed"
  out$valid <- why == ""
  out$invalid_reason <- why
  out
}

# ---- compositions ----------------------------------------------------------

COMPOSITION_SCHEMES <- c(
  sleep24  = "24 h: Sleep / Awake in bed / Out of bed (3 parts, 1440 min)",
  detail24 = "24 h: TST / SOL / WASO / EMA / Out of bed (5 parts, 1440 min)",
  tib4     = "Time in bed: TST / SOL / WASO / EMA (4 parts, closed to mean TIB)",
  custom   = "Custom: choose part columns (data already compositional)"
)

#' Multiplicative replacement of zeros (Martín-Fernández et al., 2003):
#' zeros become 0.65 x detection limit, non-zero parts shrink proportionally
#' so each row keeps its total.
replace_zeros_mult <- function(X, dl = 1) {
  X <- as.matrix(X)
  for (i in seq_len(nrow(X))) {
    z <- X[i, ] == 0
    if (!any(z, na.rm = TRUE)) next
    tot <- sum(X[i, ])
    delta <- 0.65 * dl
    X[i, z] <- delta
    X[i, !z] <- X[i, !z] * (1 - sum(z) * delta / tot)
  }
  X
}

#' Analysis data set with composition parts as columns.
#'
#' @param diary      output of derive_sleep_diary() (or any long table for "custom")
#' @param scheme     one of names(COMPOSITION_SCHEMES)
#' @param naps_as_sleep add nap minutes to sleep and remove them from out of bed
#' @param custom_parts, custom_total for scheme = "custom"
#' @param zeros      "multiplicative" or "drop"
#' @param dl         detection limit in minutes for the multiplicative replacement
#' @param keep       further columns to carry (outcomes, covariates)
#' @param shift      named integer vector: variable -> rows to shift within
#'                   person (+1 = take the next night's value, -1 = previous)
build_composition <- function(diary, scheme = "sleep24", naps_as_sleep = FALSE,
                              custom_parts = NULL, custom_total = 1440,
                              zeros = c("multiplicative", "drop"), dl = 1,
                              keep = character(0), shift = integer(0)) {
  zeros <- match.arg(zeros)
  d <- diary
  if (!"night" %in% names(d)) d$night <- stats::ave(seq_len(nrow(d)), d$ID, FUN = seq_along)
  d <- d[order(d$ID, d$night), , drop = FALSE]
  for (v in names(shift)) {
    k <- shift[[v]]
    if (!v %in% names(d) || k == 0) next
    d[[v]] <- stats::ave(d[[v]], d$ID, FUN = function(x) {
      if (k > 0) c(x[-seq_len(k)], rep(NA, k)) else c(rep(NA, -k), x[seq_len(length(x) + k)])
    })
  }
  if (scheme != "custom") d <- d[d$valid %in% TRUE, , drop = FALSE]
  nap <- if (isTRUE(naps_as_sleep) && "NAP" %in% names(d)) d$NAP else 0
  if (scheme == "sleep24") {
    X <- cbind(Sleep = d$TST + nap, WakeInBed = d$TIB - d$TST,
               OutOfBed = 1440 - d$TIB - nap)
    total <- 1440
  } else if (scheme == "detail24") {
    X <- cbind(TST = d$TST + nap, SOL = d$SOL, WASO = d$WASO,
               EMA = d$TIB - d$TST - d$SOL - d$WASO,
               OutOfBed = 1440 - d$TIB - nap)
    total <- 1440
  } else if (scheme == "tib4") {
    X <- cbind(TST = d$TST, SOL = d$SOL, WASO = d$WASO,
               EMA = d$TIB - d$TST - d$SOL - d$WASO)
    total <- round(mean(d$TIB, na.rm = TRUE))
  } else {
    if (length(custom_parts) < 2) stop("Choose at least two part columns.")
    X <- as.matrix(d[, custom_parts, drop = FALSE])
    total <- custom_total
  }
  X[abs(X) < 1e-9] <- 0
  bad <- !stats::complete.cases(X) | apply(X < 0, 1, any, na.rm = TRUE)
  n_missing <- sum(bad)
  X <- X[!bad, , drop = FALSE]
  d <- d[!bad, , drop = FALSE]
  has_zero <- apply(X == 0, 1, any)
  n_zero <- sum(has_zero)
  if (zeros == "drop") {
    X <- X[!has_zero, , drop = FALSE]
    d <- d[!has_zero, , drop = FALSE]
  } else if (n_zero) {
    X <- replace_zeros_mult(X, dl = dl)
  }
  parts <- colnames(X)
  base <- intersect(c("ID", "night", "date"), names(d))
  carry <- setdiff(intersect(keep, names(d)), c(base, parts))
  out <- cbind(d[, base, drop = FALSE], as.data.frame(X), d[, carry, drop = FALSE])
  rownames(out) <- NULL
  attr(out, "parts") <- parts
  attr(out, "total") <- total
  attr(out, "log") <- c(
    sprintf("Composition '%s': %d parts (%s), total %s min",
            scheme, length(parts), paste(parts, collapse = ", "), total),
    sprintf("Rows dropped for missing or negative parts: %d", n_missing),
    if (zeros == "drop") sprintf("Rows dropped because a part was 0: %d", n_zero)
    else sprintf("Rows with zeros replaced multiplicatively (0.65 x %s min): %d", dl, n_zero),
    sprintf("Analysis rows: %d nights from %d participant(s)",
            nrow(out), length(unique(out$ID))))
  out
}

# Closed geometric mean (compositional centre) of each part, in minutes.
comp_center <- function(X, total) {
  g <- exp(colMeans(log(as.matrix(X))))
  total * g / sum(g)
}

# ---- figures ---------------------------------------------------------------

pastel_palette <- function(n) {
  base <- c(RColorBrewer::brewer.pal(9, "Pastel1"), RColorBrewer::brewer.pal(8, "Pastel2"))
  if (n <= length(base)) base[seq_len(n)] else grDevices::colorRampPalette(base)(n)
}

darken <- function(col, f = 0.65) {
  m <- grDevices::col2rgb(col) * f
  grDevices::rgb(m[1, ], m[2, ], m[3, ], maxColorValue = 255)
}

#' Stacked bars of each person's compositional mean (minutes per day).
plot_composition_by_id <- function(ana, parts, total, colours = NULL,
                                   base_size = 12) {
  ids <- unique(ana$ID)
  rows <- lapply(ids, function(i) {
    cm <- comp_center(ana[ana$ID == i, parts, drop = FALSE], total)
    data.frame(ID = i, part = parts, minutes = as.numeric(cm))
  })
  long <- do.call(rbind, rows)
  long$part <- factor(long$part, levels = rev(parts))
  colours <- colours %||% stats::setNames(pastel_palette(length(parts)), parts)
  ggplot2::ggplot(long, ggplot2::aes(x = ID, y = minutes, fill = part)) +
    ggplot2::geom_col(colour = "grey40", linewidth = 0.2) +
    ggplot2::scale_fill_manual(values = colours, breaks = parts) +
    ggplot2::labs(x = NULL, y = "Minutes (compositional mean)", fill = NULL) +
    ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
                   legend.position = "bottom")
}

#' Night-by-night composition for one participant.
plot_composition_nights <- function(ana, parts, id, colours = NULL,
                                    base_size = 12) {
  d <- ana[ana$ID == id, , drop = FALSE]
  long <- data.frame(night = rep(d$night, length(parts)),
                     part = rep(parts, each = nrow(d)),
                     minutes = unlist(d[, parts], use.names = FALSE))
  long$part <- factor(long$part, levels = rev(parts))
  colours <- colours %||% stats::setNames(pastel_palette(length(parts)), parts)
  ggplot2::ggplot(long, ggplot2::aes(x = night, y = minutes, fill = part)) +
    ggplot2::geom_col(colour = "grey40", linewidth = 0.2, width = 0.9) +
    ggplot2::scale_fill_manual(values = colours, breaks = parts) +
    ggplot2::labs(x = "Night", y = "Minutes", fill = NULL, title = id) +
    ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(legend.position = "bottom")
}

#' Substitution estimates (summary of a multilevelcoda substitution object).
plot_substitution <- function(sub, colours = NULL, base_size = 12) {
  d <- sub$delta
  tab <- summary(sub, delta = sort(unique(c(-abs(d), abs(d)))), digits = "asis")
  tab <- as.data.frame(tab)
  parts <- sub$parts
  colours <- colours %||% stats::setNames(pastel_palette(length(parts)), parts)
  tab$From <- factor(tab$From, levels = parts)
  tab$To   <- factor(paste("to", tab$To), levels = paste("to", parts))
  p <- ggplot2::ggplot(tab, ggplot2::aes(x = Delta, y = Estimate,
                                         colour = From, fill = From)) +
    ggplot2::geom_hline(yintercept = 0, linetype = 2, linewidth = 0.3) +
    ggplot2::geom_vline(xintercept = 0, linetype = 2, linewidth = 0.3)
  if (length(unique(tab$Delta)) > 2) {
    p <- p + ggplot2::geom_ribbon(ggplot2::aes(ymin = CI_low, ymax = CI_high),
                                  alpha = 0.45, colour = NA) +
      ggplot2::geom_line(linewidth = 0.8)
  } else {
    p <- p + ggplot2::geom_pointrange(ggplot2::aes(ymin = CI_low, ymax = CI_high),
                                      position = ggplot2::position_dodge(width = 4))
  }
  p + ggplot2::scale_fill_manual(values = colours, drop = FALSE) +
    ggplot2::scale_colour_manual(values = stats::setNames(darken(colours), names(colours)),
                                 drop = FALSE) +
    ggplot2::facet_grid(Level ~ To) +
    ggplot2::labs(x = "Minutes moved from 'From' to the facet's part",
                  y = "Difference in outcome (posterior mean, 95% CrI)",
                  colour = "From", fill = "From") +
    ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(strip.background = ggplot2::element_rect(fill = "grey97"),
                   legend.position = "bottom")
}
