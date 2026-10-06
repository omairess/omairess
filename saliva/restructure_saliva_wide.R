## =============================================================================
##  restructure_saliva_wide.R
## -----------------------------------------------------------------------------
##  Salivary cortisol and melatonin (Halley 2015, Halley 2016, DDU 2017-2018,
##  DDU 2022): long format (one row per saliva sample) -> mixed format.
##
##    rows    : one per participant x sampling session (Campaign x ID x Code).
##              Between-subject factors (Location, Campaign, Sex) and
##              within-subject factors (Season, Code; Antarctic = austral
##              seasons) stay in LONG format.
##    columns : one column per clock-time bin and hormone (time of day is
##              WIDE): Cortisol_00h00 ... Cortisol_23h30,
##              Melatonin_00h00 ... Melatonin_23h30 (30-min bins by default)
##
##  Bins        : bin_min sets the width (30 = half-hour bins, 60 = hourly).
##                bin_rule "nearest" (default): column hh:mm holds the samples
##                taken within +-bin_min/2 of hh:mm (06h30 = 06:15-06:44), so
##                a sample taken a few minutes before or after a scheduled
##                time stays with that time point. Midnight (24:00 = 00:00) is
##                column 00h00 (23:45-00:14).
##                bin_rule "floor": column hh:mm holds hh:mm to hh:mm+bin_min.
##  Aggregation : several samples of one session in one bin (e.g. the 06:30
##                samples of day 1 and day 2) -> arithmetic mean. Sheet
##                "N_per_bin" gives the number of samples behind every cell.
##
##  Output sheets
##    Wide_30min  : analysis file (first sheet, one header row -> SPSS/R ready)
##    N_per_bin   : number of samples averaged in every bin
##    Long_Binned : every original sample with its clock time and bin
##    Data_Issues : entries that need checking, and how they were handled
##    Codebook    : variables and processing rules
##
##  Usage    : Rscript restructure_saliva_wide.R [input.xlsx] [output.xlsx]
##             or set in_file / out_file below and source() the script.
##  Packages : readxl, openxlsx (everything else is base R)
## =============================================================================

suppressPackageStartupMessages({
  library(readxl)
  library(openxlsx)
})

## ---- 0. Settings -------------------------------------------------------------------
bin_min  <- 30          # bin width in minutes; must divide 1440 (60 = hourly bins)
bin_rule <- "nearest"   # "nearest" or "floor" (see header)

args     <- commandArgs(trailingOnly = TRUE)
in_file  <- if (length(args) >= 1) args[1] else "Merged_Saliva_SPSSReady_1.xlsx"
out_file <- if (length(args) >= 2) args[2] else sprintf("Merged_Saliva_Wide_%dminBins.xlsx", bin_min)
in_sheet <- "Sheet1"

if (!(bin_min >= 1 && bin_min == round(bin_min) && 1440 %% bin_min == 0))
  stop("bin_min must be a whole number of minutes that divides 1440 (e.g. 15, 30, 60)")
if (!bin_rule %in% c("nearest", "floor")) stop("bin_rule must be \"nearest\" or \"floor\"")

n_bins     <- 1440 %/% bin_min
bin_start  <- (seq_len(n_bins) - 1) * bin_min                # minutes after midnight
bin_labels <- sprintf("%02dh%02d", bin_start %/% 60, bin_start %% 60)
bin_shift  <- if (bin_rule == "nearest") bin_min / 2 else 0
bin_of     <- function(s) floor((s / 60 + bin_shift) / bin_min) %% n_bins   # s: seconds after midnight
wide_sheet <- sprintf("Wide_%dmin", bin_min)

## ---- 1. Read ---------------------------------------------------------------------
## Time is read as the raw Excel serial (fraction of a day) so that cells typed
## as 1, 24 or 999 can be recognised. readxl then warns once per time cell
## ("Expecting numeric ... got a date"); only those warnings are muffled.
col_types <- c(Campaign = "text", CampaignNumber = "numeric", Order = "text",
               ID = "text", Sex_1Male_2Female = "numeric", Code = "text",
               Season = "text", Season_Summer1_Winter2 = "numeric",
               Date = "text", Time = "numeric", MinutesOfDaylight = "numeric",
               HoursOfDaylight = "numeric", Cortisol = "numeric",
               Melatonin = "numeric", BadData = "numeric", BadDataReason = "text")

hdr <- names(read_excel(in_file, sheet = in_sheet, n_max = 0))
if (!identical(hdr, names(col_types)))
  stop("Unexpected columns in ", in_file, " [", in_sheet, "]: ",
       paste(hdr, collapse = ", "))

time_col <- LETTERS[match("Time", hdr)]
raw <- withCallingHandlers(
  read_excel(in_file, sheet = in_sheet, col_types = unname(col_types)),
  warning = function(w) {
    if (grepl(paste0("^Expecting numeric in ", time_col, "[0-9]+ .*got a date"),
              conditionMessage(w)))
      invokeRestart("muffleWarning")
  })
raw <- as.data.frame(raw, stringsAsFactors = FALSE)
raw$SourceRow <- seq_len(nrow(raw)) + 1L          # row number in the Excel sheet

## ---- 2. Dates --------------------------------------------------------------------
## Date cells come in as text: Excel serials ("42388") become dates; text such as
## "D+1" (next-morning sample of a session with unknown date) goes to DateNote.
is_serial    <- !is.na(raw$Date) & grepl("^[0-9]+(\\.[0-9]+)?$", raw$Date)
raw$DateNote <- ifelse(!is.na(raw$Date) & !is_serial, raw$Date, NA_character_)
raw$Date     <- as.Date(ifelse(is_serial, floor(suppressWarnings(as.numeric(raw$Date))),
                               NA_real_), origin = "1899-12-30")

## ---- 3. Clock time and time bin --------------------------------------------------
## A valid clock time is an Excel serial in [0, 1). Anything else is handled as:
##   1               -> 24:00 typed as a time         -> 00:00
##   integer 2 to 24 -> hours typed without a colon   -> (value mod 24):00
##   other values    -> not a clock time (e.g. 999 = missing-value code) -> no bin
tv          <- raw$Time
sec         <- rep(NA_real_, length(tv))
t_problem   <- rep(NA_character_, length(tv))
t_handling  <- rep(NA_character_, length(tv))
col_of      <- function(s) paste("column", bin_labels[bin_of(s) + 1])

is_clock <- !is.na(tv) & tv >= 0 & tv < 1
sec[is_clock] <- round(tv[is_clock] * 86400)

is_one <- !is.na(tv) & tv == 1
sec[is_one]        <- 0
t_problem[is_one]  <- "Time cell holds 1 (Excel serial for 24:00)"
t_handling[is_one] <- paste("Read as 24:00 = 00:00 (midnight) ->", col_of(0))

is_hours <- !is.na(tv) & tv > 1 & tv <= 24 & tv == round(tv)
sec[is_hours]        <- (tv[is_hours] %% 24) * 3600
t_problem[is_hours]  <- sprintf("Time cell holds %g (hours typed without a colon)", tv[is_hours])
t_handling[is_hours] <- sprintf("Read as %g:00 = %02d:00 -> %s", tv[is_hours],
                                as.integer(tv[is_hours] %% 24), col_of(sec[is_hours]))

is_invalid <- !is.na(tv) & !(is_clock | is_one | is_hours)
t_problem[is_invalid]  <- sprintf("Time cell holds %g: not a clock time (probably a missing-value code)",
                                  tv[is_invalid])
t_handling[is_invalid] <- "Sample not assigned to a bin (values kept in Long_Binned)"

t_problem[is.na(tv)]  <- "Time missing"
t_handling[is.na(tv)] <- "Sample not assigned to a bin (values kept in Long_Binned)"

sec             <- sec %% 86400                    # 24:00:00 after rounding -> 00:00:00
raw$Time        <- sec / 86400                     # cleaned Excel time (fraction of a day)
raw$ClockTime_h <- sec / 3600
raw$Bin         <- bin_of(sec)                     # 0-based bin index, NA without time
raw$BinLabel    <- bin_labels[raw$Bin + 1]
raw$BinTime_h   <- bin_start[raw$Bin + 1] / 60
raw$TimeNote    <- ifelse(is.na(t_problem), NA_character_,
                          paste0(t_problem, ". ", t_handling))

## ---- 4. Derived identifiers ------------------------------------------------------
raw$Location <- sub("[0-9]+$", "", raw$Campaign)   # Halley2016 -> Halley, DDU2022 -> DDU
if (!all(raw$Location %in% c("Halley", "DDU")))
  stop("Unknown campaign prefix: ",
       paste(unique(raw$Campaign[!raw$Location %in% c("Halley", "DDU")]), collapse = ", "))
raw$Location_Halley1_DDU2 <- ifelse(raw$Location == "Halley", 1, 2)
raw$SessionID <- paste(raw$Campaign, raw$ID, raw$Code, sep = "_")

raw$InBins_Cortisol  <- as.integer(!is.na(raw$Cortisol)  & !is.na(raw$Bin))
raw$InBins_Melatonin <- as.integer(!is.na(raw$Melatonin) & !is.na(raw$Bin))

## ---- 5. One row per participant x session ----------------------------------------
## SessionDate = earliest recorded date in the session (day 1). Daylight is taken
## from the rows recorded on that date (or, without any date, from the first
## non-missing value of the session).
first_non_na <- function(x) { x <- x[!is.na(x)]; if (length(x)) x[1] else x[NA_integer_] }

constant_vars <- c("Campaign", "CampaignNumber", "Location", "Location_Halley1_DDU2",
                   "ID", "Sex_1Male_2Female", "Code", "Season", "Season_Summer1_Winter2")

rows_by_session <- split(seq_len(nrow(raw)),
                         factor(raw$SessionID, levels = unique(raw$SessionID)))

sess <- do.call(rbind, lapply(rows_by_session, function(i) {
  for (v in constant_vars)
    if (length(unique(raw[[v]][i])) > 1)
      stop(v, " is not constant within session ", raw$SessionID[i[1]])
  d     <- raw$Date[i]
  d0    <- if (all(is.na(d))) as.Date(NA) else min(d, na.rm = TRUE)
  on_d0 <- if (is.na(d0)) i else i[!is.na(d) & d == d0]
  why   <- unique(stats::na.omit(raw$BadDataReason[i]))
  data.frame(SessionID         = raw$SessionID[i[1]],
             raw[i[1], constant_vars],
             SessionDate       = d0,
             MinutesOfDaylight = first_non_na(raw$MinutesOfDaylight[on_d0]),
             HoursOfDaylight   = first_non_na(raw$HoursOfDaylight[on_d0]),
             BadData           = if (any(raw$BadData[i] == 1, na.rm = TRUE)) 1 else NA_real_,
             BadDataReason     = if (length(why)) paste(why, collapse = "; ") else NA_character_,
             N_Samples         = length(i),
             N_Cortisol        = sum(raw$InBins_Cortisol[i]),
             N_Melatonin       = sum(raw$InBins_Melatonin[i]),
             stringsAsFactors  = FALSE, row.names = NULL)
}))

## Order: campaign, ID (natural order: SH1, SH2, ..., SH11), session date
## (sessions with unknown date last), Code.
id_stem <- sub("[0-9]+$", "", sess$ID)
id_num  <- suppressWarnings(as.numeric(sub("^[^0-9]*", "", sess$ID)))
sess    <- sess[order(sess$CampaignNumber, id_stem, id_num, sess$SessionDate, sess$Code,
                      na.last = TRUE), ]
rownames(sess) <- NULL

## ---- 6. Time bins: means and counts ----------------------------------------------
bin_table <- function(value, FUN) {
  use <- !is.na(value) & !is.na(raw$Bin)
  m <- tapply(value[use],
              list(factor(raw$SessionID[use], levels = sess$SessionID),
                   factor(raw$Bin[use],       levels = seq_len(n_bins) - 1)),
              FUN)
  matrix(as.numeric(m), nrow = nrow(sess))        # rows follow sess, columns follow bins
}
cort_mean <- bin_table(raw$Cortisol,  mean)
mel_mean  <- bin_table(raw$Melatonin, mean)
cort_n    <- bin_table(raw$Cortisol,  length); cort_n[is.na(cort_n)] <- 0
mel_n     <- bin_table(raw$Melatonin, length); mel_n[is.na(mel_n)]   <- 0
colnames(cort_mean) <- paste0("Cortisol_",    bin_labels)
colnames(mel_mean)  <- paste0("Melatonin_",   bin_labels)
colnames(cort_n)    <- paste0("N_Cortisol_",  bin_labels)
colnames(mel_n)     <- paste0("N_Melatonin_", bin_labels)

## ---- 7. Checks -> Data_Issues ----------------------------------------------------
fmt_hm <- function(s) ifelse(is.na(s), "unknown",
                             sprintf("%02d:%02d", (s %/% 3600) %% 24, (s %% 3600) %/% 60))
fmt_rows <- function(r) {
  r <- sort(unique(r)); runs <- split(r, cumsum(c(1, diff(r) != 1)))
  paste(vapply(runs, function(x) if (length(x) > 1) paste0(x[1], "-", x[length(x)]) else
    as.character(x), character(1)), collapse = ", ")
}
issue_rows <- list()
add_issue <- function(check, session, rows, problem, handling)
  issue_rows[[length(issue_rows) + 1]] <<- data.frame(
    Check = check, SessionID = session, SourceRows = rows,
    Problem = problem, Handling = handling, stringsAsFactors = FALSE)

## (a) Time cells that are not a clock time
for (r in which(!is.na(t_problem)))
  add_issue("Time value", raw$SessionID[r], as.character(raw$SourceRow[r]),
            sprintf("Order %s: %s", ifelse(is.na(raw$Order[r]), "-", raw$Order[r]), t_problem[r]),
            t_handling[r])

## (b) Recorded (austral) Season vs the daylight recorded at the station on the
##     sampling date: Summer with < 12 h or Winter with > 12 h of daylight
mismatch <- which(!is.na(raw$MinutesOfDaylight) &
                  ((raw$Season == "Summer" & raw$MinutesOfDaylight < 720) |
                   (raw$Season == "Winter" & raw$MinutesOfDaylight > 720)))
for (s in unique(raw$SessionID[mismatch])) {
  r <- mismatch[raw$SessionID[mismatch] == s]
  summer_label <- raw$Season[r[1]] == "Summer"
  add_issue("Season vs photoperiod", s, fmt_rows(raw$SourceRow[r]),
            sprintf("Code = %s, Season = %s, but %d of %d samples are dated %s, with %s min of daylight at the station (%s 12 h: austral %s half-year)",
                    raw$Code[r[1]], raw$Season[r[1]], length(r), sum(raw$SessionID == s),
                    paste(unique(format(raw$Date[r])), collapse = ", "),
                    paste(unique(raw$MinutesOfDaylight[r]), collapse = ", "),
                    ifelse(summer_label, "<", ">"), ifelse(summer_label, "winter", "summer")),
            "Not changed: check whether Code/Season or Date is wrong. The session is filed under the recorded Season.")
}

## (c) Recorded dates of one session more than 1 day apart
for (s in sess$SessionID) {
  i <- which(raw$SessionID == s)
  d <- raw$Date[i]
  if (sum(!is.na(d)) > 1 && as.numeric(diff(range(d, na.rm = TRUE))) > 1)
    add_issue("Date span", s, fmt_rows(raw$SourceRow[i]),
              sprintf("Recorded dates range from %s to %s (%d days)",
                      min(d, na.rm = TRUE), max(d, na.rm = TRUE),
                      as.integer(diff(range(d, na.rm = TRUE)))),
              sprintf("SessionDate = earliest date (%s); daylight taken from the rows on that date",
                      min(d, na.rm = TRUE)))
}

## (d) Two samples of one session with identical date and clock time
day_key <- ifelse(is.na(raw$Date), ifelse(is.na(raw$DateNote), "day?", raw$DateNote),
                  format(raw$Date))
stamp   <- ifelse(is.na(sec), NA, paste(raw$SessionID, day_key, sec))
for (st in unique(stamp[!is.na(stamp) & duplicated(stamp)])) {
  r <- which(stamp == st)
  add_issue("Duplicate time", raw$SessionID[r[1]], fmt_rows(raw$SourceRow[r]),
            sprintf("Orders %s have the same recorded date and time (%s %s)",
                    paste(raw$Order[r], collapse = " and "), day_key[r[1]], fmt_hm(sec[r[1]])),
            sprintf(paste("Binned as recorded (column %s, mean of %d samples). If one is a",
                          "next-morning sample whose date was not advanced, binning is",
                          "unaffected; otherwise check the time entries"),
                    raw$BinLabel[r[1]], length(r)))
}

## (e) DDU protocols (sheet "Timing" of the input file): fixed-clock samples
##     T6 = 12:00, T9 = 15:00, T12 = 18:00, T15 = 21:00, T16 = 22:00,
##     T17 = 23:00, T18 = 24:00. Flag recorded times >= 90 min away.
sched_h  <- c("6" = 12, "9" = 15, "12" = 18, "15" = 21, "16" = 22, "17" = 23, "18" = 0)
code_num <- ifelse(raw$Location == "DDU" & grepl("^T?[0-9]+$", raw$Order),
                   sub("^T", "", raw$Order), NA)
fixed    <- which(!is.na(code_num) & code_num %in% names(sched_h) & !is.na(sec))
dev_min  <- ((sec[fixed] / 60 - sched_h[code_num[fixed]] * 60 + 720) %% 1440) - 720
for (k in which(abs(dev_min) >= 90)) {
  r <- fixed[k]
  add_issue("Protocol time", raw$SessionID[r], as.character(raw$SourceRow[r]),
            sprintf("Order %s (protocol %02d:00) recorded at %s (%+d min)",
                    raw$Order[r], as.integer(sched_h[code_num[r]]), fmt_hm(sec[r]),
                    as.integer(round(dev_min[k]))),
            sprintf("Binned as recorded (column %s); check the time entry", raw$BinLabel[r]))
}

## (f) DDU protocols: wake sample (T0) later than the next sample (T1) on the same day
for (s in unique(raw$SessionID[raw$Location == "DDU"])) {
  r0 <- which(raw$SessionID == s & code_num == "0" & !is.na(sec))
  r1 <- which(raw$SessionID == s & code_num == "1" & !is.na(sec))
  if (length(r0) == 1 && length(r1) == 1 && identical(day_key[r0], day_key[r1]) && sec[r0] > sec[r1])
    add_issue("Protocol order", s, as.character(raw$SourceRow[r0]),
              sprintf("Wake sample %s recorded at %s, later than %s at %s on the same day",
                      raw$Order[r0], fmt_hm(sec[r0]), raw$Order[r1], fmt_hm(sec[r1])),
              sprintf("Binned as recorded (column %s); check whether the time or the Order label is wrong",
                      raw$BinLabel[r0]))
}

issues <- if (length(issue_rows)) do.call(rbind, issue_rows) else
  data.frame(Check = character(), SessionID = character(), SourceRows = character(),
             Problem = character(), Handling = character())
check_order <- c("Time value", "Season vs photoperiod", "Date span", "Duplicate time",
                 "Protocol time", "Protocol order")
issues <- issues[order(match(issues$Check, check_order), match(issues$SessionID, sess$SessionID)), ]
issues <- data.frame(IssueID = sprintf("I%02d", seq_len(nrow(issues))), issues,
                     stringsAsFactors = FALSE, row.names = NULL)

sess$DataIssues <- vapply(sess$SessionID, function(s) {
  id <- issues$IssueID[issues$SessionID == s]
  if (length(id)) paste(id, collapse = "; ") else NA_character_
}, character(1), USE.NAMES = FALSE)

## ---- 8. Assemble sheets ----------------------------------------------------------
wide   <- data.frame(sess, cort_mean, mel_mean, check.names = FALSE, stringsAsFactors = FALSE)
counts <- data.frame(sess[, c("SessionID", "Campaign", "ID", "Code", "Season")],
                     cort_n, mel_n, check.names = FALSE, stringsAsFactors = FALSE)
long   <- raw[, c("SourceRow", "SessionID", "Campaign", "CampaignNumber", "Location",
                  "Location_Halley1_DDU2", "ID", "Sex_1Male_2Female", "Code", "Season",
                  "Season_Summer1_Winter2", "Order", "Date", "DateNote", "Time",
                  "ClockTime_h", "BinLabel", "BinTime_h", "TimeNote", "MinutesOfDaylight",
                  "HoursOfDaylight", "Cortisol", "Melatonin", "BadData", "BadDataReason",
                  "InBins_Cortisol", "InBins_Melatonin")]

## SPSS variable-name rules: letter first, letters/digits/underscore, <= 64 bytes
for (nm in list(names(wide), names(counts), names(long)))
  stopifnot(all(grepl("^[A-Za-z][A-Za-z0-9_]*$", nm)), all(nchar(nm) <= 64), !anyDuplicated(nm))

## facts for the codebook, computed from the data (stay correct after source fixes)
camp_map  <- unique(raw[order(raw$CampaignNumber), c("Campaign", "CampaignNumber")])
id_camps  <- tapply(raw$Campaign, raw$ID, function(x) length(unique(x)))
dl_varies <- names(which(tapply(raw$MinutesOfDaylight, raw$SessionID,
                                function(x) length(unique(stats::na.omit(x)))) > 1))
reasons   <- unique(stats::na.omit(raw$BadDataReason))
code_txt  <- vapply(split(raw, raw$CampaignNumber), function(d) {
  per_id <- tapply(d$Code, d$ID, function(x) length(unique(x)))
  sprintf("%s: %s (%d of %d IDs with more than one session)", d$Campaign[1],
          paste(sort(unique(d$Code)), collapse = " / "), sum(per_id > 1), length(per_id))
}, character(1))

fmt_min  <- function(m) sprintf("%02d:%02d", (m %% 1440) %/% 60, (m %% 1440) %% 60)
span_txt <- function(mark) {                          # clock times that go to column 'mark'
  lo <- ceiling(mark - bin_shift); hi <- ceiling(mark - bin_shift + bin_min) - 1
  paste0(fmt_min(lo), "-", fmt_min(hi))
}
ex_mark  <- 6 * 60 + 30 - (6 * 60 + 30) %% bin_min    # example column near 06:30
bin_txt  <- sprintf(paste("%d columns of %d min per hormone, labelled hh'h'mm. Column %s holds the samples",
                          "recorded %s%s. Midnight (24:00 = 00:00) is column %s (%s). Samples taken after",
                          "midnight (bedtime) sit in the first columns, as on a 0-24 h time-of-day axis."),
                    n_bins, bin_min, bin_labels[ex_mark / bin_min + 1], span_txt(ex_mark),
                    if (bin_rule == "nearest")
                      ", i.e. within half a bin of that time, so a sample taken a few minutes before or after a scheduled time stays with that time point"
                    else ", i.e. from that time to the start of the next column",
                    bin_labels[1], span_txt(0))
bins_first_last <- function(prefix) paste0(prefix, bin_labels[1], " ... ", prefix, bin_labels[n_bins])

codebook <- rbind(
  c("File", "Source", sprintf("%s, sheet '%s' (%d samples, %d columns)", basename(in_file),
                              in_sheet, nrow(raw), length(hdr))),
  c("File", "Created", sprintf("%s by restructure_saliva_wide.R (bin_min = %d, bin_rule = '%s'; R %s, readxl %s, openxlsx %s)",
                               format(Sys.Date()), bin_min, bin_rule, getRversion(),
                               packageVersion("readxl"), packageVersion("openxlsx"))),
  c("Structure", sprintf("Row unit (%s)", wide_sheet), sprintf("One row per participant x sampling session (Campaign x ID x Code): %d rows. Between-subject (Location, Campaign, Sex) and within-subject (Season, Code) variables are in long format; time of day is in wide format.", nrow(wide))),
  c("Structure", "Time bins", bin_txt),
  c("Structure", "Several samples in one bin", "Cell = arithmetic mean of all samples of that session in that bin, regardless of day (e.g. the 06:30 samples of day 1 and day 2 in DDU 2017-2018). Sheet N_per_bin gives the number of samples behind every cell (highlighted when >= 2)."),
  c("Structure", "Empty cell", "No sample (or no assay value) in that bin for that session. Read as missing by SPSS (system-missing) and R (NA)."),
  c("Rules", "SessionDate", "Earliest recorded date of the session (day 1). Empty when no date was recorded."),
  c("Rules", "MinutesOfDaylight / HoursOfDaylight", paste0("Values recorded on SessionDate; sessions without a date: first non-missing value in the session. ",
    if (length(dl_varies)) paste0("Daylight differs within: ", paste(dl_varies, collapse = ", "), " (see Data_Issues).") else
      "Daylight is constant within every session.")),
  c("Rules", "BadData / BadDataReason", paste0("1 / reason if any sample of the session was flagged in the source file; empty otherwise. Hormone values of flagged samples are kept. Reasons in the source: ",
    if (length(reasons)) paste(reasons, collapse = "; ") else "none", ".")),
  c("Rules", "Time cells that are not clock times", "Excel value 1 (= 24:00) is read as 00:00; whole numbers 2-24 (hours typed without a colon) as that hour (24 -> 00:00); any other value (e.g. 999) as a missing time, and the sample is not binned. Each case is listed in Data_Issues."),
  c("Rules", "Row order", "CampaignNumber, ID (natural order), SessionDate (unknown dates last), Code."),
  c("Rules", "Data values", sprintf("Nothing in the source was corrected. Sessions whose labels or times look inconsistent are listed in Data_Issues and in the DataIssues column of %s.", wide_sheet)),
  c("Variables", "SessionID", "Campaign_ID_Code; unique key of a session (row)."),
  c("Variables", "Campaign / CampaignNumber", paste0("As in source: ", paste(camp_map$Campaign, "=", camp_map$CampaignNumber, collapse = ", "), ". ",
    if (any(id_camps > 1)) paste("IDs in more than one campaign:", paste(names(id_camps)[id_camps > 1], collapse = ", ")) else
      "Between-subject: no ID occurs in two campaigns.")),
  c("Variables", "Location / Location_Halley1_DDU2", "Derived from Campaign: Halley (Halley VI) = 1, DDU (Dumont d'Urville) = 2. Between-subject."),
  c("Variables", "ID", "Participant code, as in source."),
  c("Variables", "Sex_1Male_2Female", "As in source (constant within ID)."),
  c("Variables", "Code", paste0("Antarctic season of the participant's stay, as in source. ", paste(code_txt, collapse = "; "), ".")),
  c("Variables", "Season / Season_Summer1_Winter2", "Antarctic (austral) season, as in source (Summer = 1, Winter = 2)."),
  c("Variables", "N_Samples", "Number of rows (saliva samples) of the session in the source file."),
  c("Variables", "N_Cortisol / N_Melatonin", "Number of non-missing values of the session that were assigned to a bin."),
  c("Variables", "DataIssues", "IssueIDs (sheet Data_Issues) that concern the session; empty when none."),
  c("Variables", bins_first_last("Cortisol_"), "Mean salivary cortisol in that time bin, same units as the source file."),
  c("Variables", bins_first_last("Melatonin_"), "Mean salivary melatonin in that time bin, same units as the source file."),
  c("Other sheets", "N_per_bin", sprintf("Same rows as %s: number of cortisol and melatonin samples behind every bin (0 = empty cell).", wide_sheet)),
  c("Other sheets", "Long_Binned", sprintf("All source rows in source order (SourceRow = Excel row in the source sheet) with Time (cleaned), ClockTime_h (decimal hours), BinLabel (column suffix in %s), BinTime_h (time of that column in decimal hours), TimeNote, and InBins_* = 1 when the value entered %s. Date text such as 'D+1' (next morning, date unknown) is in DateNote.", wide_sheet, wide_sheet)),
  c("Other sheets", "Data_Issues", "Rule-based checks: time cells that are not clock times; recorded (austral) Season vs the daylight recorded at the station on the sampling date (Summer with < 12 h or Winter with > 12 h); recorded dates > 1 day apart within a session; duplicate date+time within a session; DDU fixed-clock protocol samples (T6-T18, sheet Timing) recorded >= 90 min off; DDU wake sample later than the next sample.")
)
codebook <- data.frame(Section = codebook[, 1], Item = codebook[, 2], Description = codebook[, 3],
                       stringsAsFactors = FALSE)

## ---- 9. Write workbook -----------------------------------------------------------
options(openxlsx.dateFormat = "yyyy-mm-dd")
wb <- createWorkbook(creator = "restructure_saliva_wide.R")
modifyBaseFont(wb, fontSize = 10, fontColour = "#000000", fontName = "Arial")
st_head  <- createStyle(textDecoration = "bold", fgFill = "#D9E1F2", border = "bottom")
st_num3  <- createStyle(numFmt = "0.000")
st_num2  <- createStyle(numFmt = "0.00")
st_time  <- createStyle(numFmt = "hh:mm")
st_multi <- createStyle(fgFill = "#FCE4D6")
st_wrap  <- createStyle(wrapText = TRUE, valign = "top")

put_sheet <- function(name, df, freeze_col = 1) {
  addWorksheet(wb, name)
  writeData(wb, name, df, headerStyle = st_head, keepNA = FALSE)
  freezePane(wb, name, firstActiveRow = 2, firstActiveCol = freeze_col + 1)
  addFilter(wb, name, rows = 1, cols = seq_len(ncol(df)))
}
fit_widths <- function(name, df, cols = seq_len(ncol(df)), max_w = 45) {
  w <- vapply(cols, function(j) {
    v   <- df[[j]]
    txt <- if (inherits(v, "Date")) "yyyy-mm-dd" else
           if (is.numeric(v)) formatC(v, digits = 3, format = "fg") else as.character(v)
    max(nchar(names(df)[j]), nchar(txt), na.rm = TRUE)
  }, numeric(1))
  setColWidths(wb, name, cols = cols, widths = pmin(ceiling(w * 1.1) + 2, max_w))
}
style_cols <- function(name, df, cols, style)
  if (length(cols) && nrow(df))
    addStyle(wb, name, style, rows = seq_len(nrow(df)) + 1, cols = cols,
             gridExpand = TRUE, stack = TRUE)

put_sheet(wide_sheet, wide)
bin_cols <- grep("^(Cortisol|Melatonin)_[0-9]{2}h[0-9]{2}$", names(wide))
style_cols(wide_sheet, wide, bin_cols, st_num3)
style_cols(wide_sheet, wide, match("HoursOfDaylight", names(wide)), st_num2)
fit_widths(wide_sheet, wide, cols = seq_len(ncol(sess)))
setColWidths(wb, wide_sheet, cols = bin_cols, widths = 17)

put_sheet("N_per_bin", counts)
n_cols <- grep("^N_(Cortisol|Melatonin)_[0-9]{2}h[0-9]{2}$", names(counts))
conditionalFormatting(wb, "N_per_bin", cols = n_cols, rows = seq_len(nrow(counts)) + 1,
                      rule = ">=2", style = st_multi)
fit_widths("N_per_bin", counts, cols = 1:5)
setColWidths(wb, "N_per_bin", cols = n_cols, widths = 19)

put_sheet("Long_Binned", long, freeze_col = 2)
style_cols("Long_Binned", long, match("Time", names(long)), st_time)
style_cols("Long_Binned", long, match(c("ClockTime_h", "BinTime_h", "HoursOfDaylight"), names(long)), st_num2)
fit_widths("Long_Binned", long)
setColWidths(wb, "Long_Binned", cols = match("TimeNote", names(long)), widths = 60)

put_sheet("Data_Issues", issues)
style_cols("Data_Issues", issues, seq_len(ncol(issues)), st_wrap)
setColWidths(wb, "Data_Issues", cols = 1:6, widths = c(8, 20, 30, 12, 70, 60))

put_sheet("Codebook", codebook, freeze_col = 0)
style_cols("Codebook", codebook, 1:3, st_wrap)
setColWidths(wb, "Codebook", cols = 1:3, widths = c(12, 34, 110))

saveWorkbook(wb, out_file, overwrite = TRUE)

## openxlsx keeps a relationship and a content-type entry for a drawing part of
## every sheet, also when it writes no such part (no images). Excel ignores them,
## stricter readers (e.g. Python openpyxl) refuse the file. Drop the dangling
## entries; nothing else in the package is touched.
drop_dangling_parts <- function(xlsx) {
  xlsx <- normalizePath(xlsx)
  tmp  <- tempfile("xlsx_")
  dir.create(tmp)
  on.exit(unlink(tmp, recursive = TRUE), add = TRUE)
  zip::unzip(xlsx, exdir = tmp)
  read_xml_text  <- function(f) paste(readLines(f, warn = FALSE, encoding = "UTF-8"), collapse = "")
  write_xml_text <- function(x, f) writeLines(x, f, useBytes = TRUE)
  changed <- FALSE
  for (rel in list.files(file.path(tmp, "xl", "worksheets", "_rels"), full.names = TRUE)) {
    x    <- read_xml_text(rel)
    tags <- regmatches(x, gregexpr("<Relationship [^>]*/>", x))[[1]]
    tgt  <- sub('.*Target="([^"]+)".*', "\\1", tags)
    gone <- !grepl("TargetMode=\"External\"", tags) &
            !file.exists(file.path(tmp, "xl", "worksheets", tgt))
    for (t in tags[gone]) x <- sub(t, "", x, fixed = TRUE)
    if (any(gone)) { write_xml_text(x, rel); changed <- TRUE }
  }
  ct   <- file.path(tmp, "[Content_Types].xml")
  x    <- read_xml_text(ct)
  tags <- regmatches(x, gregexpr("<Override [^>]*/>", x))[[1]]
  gone <- !file.exists(file.path(tmp, sub('.*PartName="/([^"]+)".*', "\\1", tags)))
  for (t in tags[gone]) x <- sub(t, "", x, fixed = TRUE)
  if (any(gone)) { write_xml_text(x, ct); changed <- TRUE }
  if (changed) {
    parts <- list.files(tmp, recursive = TRUE, all.files = TRUE, no.. = TRUE)
    zip::zip(xlsx, files = parts, root = tmp, mode = "mirror")
  }
  invisible(changed)
}
drop_dangling_parts(out_file)

## ---- 10. Console summary ---------------------------------------------------------
cat(sprintf("\nRead  : %s (%d samples)\nWrote : %s\n", in_file, nrow(raw), out_file))
cat(sprintf("%s: %d sessions x %d columns (%d session variables + 2 x %d bins of %d min, rule '%s')\n",
            wide_sheet, nrow(wide), ncol(wide), ncol(sess), n_bins, bin_min, bin_rule))
print(table(Campaign = sess$Campaign, Code = sess$Code))
cat(sprintf("\nFilled bins  : cortisol %d of %d (%.0f%%), melatonin %d of %d (%.0f%%)\n",
            sum(cort_n > 0), length(cort_n), 100 * mean(cort_n > 0),
            sum(mel_n > 0), length(mel_n), 100 * mean(mel_n > 0)))
cat(sprintf("Bins = mean of >= 2 samples: cortisol %d, melatonin %d\n",
            sum(cort_n >= 2), sum(mel_n >= 2)))
cat(sprintf("Samples not binned (no valid time): %d\n", sum(is.na(raw$Bin))))
cat(sprintf("Data issues listed: %d\n", nrow(issues)))
print(table(issues$Check)[intersect(check_order, issues$Check)])
