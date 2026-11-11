# ==============================================================================
# server/02e_helpers_numeric.R -- curve columns as numbers, or as declared gaps
#
# WHY THIS EXISTS (P19.1)
# -----------------------
# The import step built the analysis matrix with data.matrix(). On a character
# column data.matrix() does not fail: it converts the column to a FACTOR and
# returns the factor's integer codes. So a curve column that R had read as text
# -- a decimal comma ("3,5"), a "<LOQ" flag, a stray unit -- went into every
# analysis as 1, 2, 3 ... in alphabetical order of its strings, with no error
# and no warning. The as.numeric() fallback written next to it sat in an error
# handler that data.matrix() never triggers.
#
# Two pure functions replace it:
#
#   dance_fix_decimal_commas()  repairs whole columns that are plainly numbers
#                               written with a decimal comma, and says which;
#   dance_numeric_matrix()      converts the selected curve columns with
#                               as.numeric(), never through a factor, and counts
#                               every cell that held text it could not read.
#
# Neither touches app state, so the tests exercise them directly.
# ==============================================================================

# Strings that mean "no value here". Compared after trimming, case-insensitive.
DANCE_MISSING_TOKENS <- c("", "na", "n/a", "nan", "null", "#n/a", ".", "-", "--", "?")

# TRUE where a character cell is a declared gap rather than unreadable text.
dance_is_missing_token <- function(x) {
  is.na(x) | tolower(trimws(as.character(x))) %in% DANCE_MISSING_TOKENS
}

# Columns of `df` that are numbers written with a DECIMAL COMMA, turned into
# numeric columns. A column qualifies only when every non-missing cell is an
# optionally signed number whose decimal mark, if it has one, is a comma, and at
# least one cell does carry a comma. `mode`:
#   "auto"   convert the columns that qualify (the default);
#   "comma"  the same test, stated by the user rather than detected;
#   "point"  convert nothing.
# Returns the data frame and the names of the converted columns.
dance_fix_decimal_commas <- function(df, mode = c("auto", "comma", "point")) {
  mode <- match.arg(mode)
  out <- list(data = df, converted = character(0))
  if (identical(mode, "point") || !is.data.frame(df) || !ncol(df)) return(out)
  pat <- "^[+-]?[0-9]+(,[0-9]+)?([eE][+-]?[0-9]+)?$"
  for (nm in names(df)) {
    v <- df[[nm]]
    if (is.factor(v)) v <- as.character(v)
    if (!is.character(v)) next
    s <- trimws(v)
    gap <- dance_is_missing_token(s)
    if (all(gap)) next
    if (!all(grepl(pat, s[!gap])) || !any(grepl(",", s[!gap], fixed = TRUE))) next
    num <- rep(NA_real_, length(s))
    num[!gap] <- as.numeric(sub(",", ".", s[!gap], fixed = TRUE))
    df[[nm]] <- num
    out$converted <- c(out$converted, nm)
  }
  out$data <- df
  out
}

# The analysis matrix: one row per row of `df`, one column per column of `df`,
# storage mode double. `dec` is the import's decimal mark ("auto", "comma" or
# "point"). A numeric or logical column is taken as it is. A
# character or factor column is read with as.numeric() on its STRINGS; a cell
# that is a declared gap becomes NA, and a cell holding anything else that is
# not a number becomes NA too and is COUNTED, so the caller can say so.
#
# Returns
#   matrix     the numeric matrix, column names kept
#   n_bad      named integer vector, unreadable cells per column (zeros kept)
#   examples   named character vector, up to three distinct unreadable strings
#              per column that has any
#   text_cols  the columns that arrived as text at all
dance_numeric_matrix <- function(df, dec = c("auto", "comma", "point")) {
  dec <- match.arg(dec)
  comma_pat <- "^[+-]?[0-9]+,[0-9]+([eE][+-]?[0-9]+)?$"
  df <- as.data.frame(df, stringsAsFactors = FALSE, check.names = FALSE)
  n <- nrow(df); p <- ncol(df)
  mat <- matrix(NA_real_, n, p, dimnames = list(NULL, names(df)))
  n_bad <- stats::setNames(integer(p), names(df))
  examples <- character(0); text_cols <- character(0)
  for (j in seq_len(p)) {
    v <- df[[j]]
    if (is.numeric(v) || is.logical(v)) { mat[, j] <- as.numeric(v); next }
    s <- trimws(as.character(v))
    text_cols <- c(text_cols, names(df)[j])
    gap <- dance_is_missing_token(s)
    num <- rep(NA_real_, n)
    num[!gap] <- suppressWarnings(as.numeric(s[!gap]))
    # a decimal comma in a column that also carries a flag ("<LOQ"): the whole-
    # column repair above leaves such a column as text, so the cells that are
    # plainly comma decimals are read here, unless the decimal mark is "Point"
    if (!identical(dec, "point")) {
      cm <- !gap & is.na(num) & grepl(comma_pat, s)
      num[cm] <- as.numeric(sub(",", ".", s[cm], fixed = TRUE))
    }
    bad <- !gap & is.na(num)
    n_bad[j] <- sum(bad)
    if (any(bad))
      examples[names(df)[j]] <- paste(sprintf("'%s'", utils::head(unique(s[bad]), 3)),
                                      collapse = ", ")
    mat[, j] <- num
  }
  list(matrix = mat, n_bad = n_bad, examples = examples, text_cols = text_cols)
}

# One sentence per affected column, for the notification and the console.
dance_numeric_matrix_message <- function(res) {
  bad <- res$n_bad[res$n_bad > 0]
  if (!length(bad)) return(NULL)
  sprintf("%d cell%s in %d curve column%s held text that is not a number and %s read as missing: %s.",
          sum(bad), if (sum(bad) == 1) "" else "s",
          length(bad), if (length(bad) == 1) "" else "s",
          if (sum(bad) == 1) "was" else "were",
          paste(sprintf("%s (%d: %s)", names(bad), bad, res$examples[names(bad)]),
                collapse = "; "))
}
