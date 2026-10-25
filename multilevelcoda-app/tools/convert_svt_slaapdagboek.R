# ============================================================================
# convert_svt_slaapdagboek.R
#
# Restructures an export of the SVT sleep-diary app (one row per entry, the
# night stored as text in 'Slaaptijd') into one row per participant-night
# that the CoDA sleep-diary Shiny app reads directly
# (layout "One row per night", ID from the column "Participant").
#
# 'Slaaptijd' holds:  start, <bedtime>, end, <rise time>, sleepIn, <n>,
#   values, time, 18:00, status, not_in_bed, time, 18:15, status, ...,
#   totals, sl, <min>, wans, <min>, twt, <min>, tib, <min>, tst, <min>, se, <%>
# The grid has 15-min epochs (statuses: not_in_bed, in_bed, asleep, awake,
# out_of_bed, nap) starting at 18:00 or 20:00 on 'Dag', the evening the
# night begins. The app's totals equal the grid (checked below):
#   tst = asleep epochs; tib = in-bed window (in_bed + asleep + awake +
#   out_of_bed); twt = sl + wans; tib = tst + twt. 'wans' is ALL wake after
#   sleep onset, including the wake between final awakening and getting up.
#   That last part is split off here as EMA, so WASO = wans - EMA, which is
#   the Consensus Sleep Diary definition (Carney et al., 2012, Sleep 35:287).
#
# Usage:  Rscript convert_svt_slaapdagboek.R <export.xlsx> <output.xlsx>
# ============================================================================

suppressPackageStartupMessages({
  library(readxl)
  library(writexl)
})

args <- commandArgs(trailingOnly = TRUE)
in_file  <- if (length(args) >= 1) args[1] else "SVT slaapdagboeken toestemming.xlsx"
out_file <- if (length(args) >= 2) args[2] else "SVT_slaapdagboek_per_nacht.xlsx"

raw <- read_excel(in_file, sheet = 1, col_types = "text")
names(raw) <- trimws(names(raw))
stopifnot(all(c("Owner Id", "Slaaptijd", "Dag") %in% names(raw)))

# ---- parse one 'Slaaptijd' cell ---------------------------------------------
to_min <- function(hhmm) {
  p <- strsplit(hhmm, ":", fixed = TRUE)
  vapply(p, function(z) as.numeric(z[1]) * 60 + as.numeric(z[2]), numeric(1))
}
fmt_clock <- function(m) ifelse(is.na(m), NA_character_,
  sprintf("%02d:%02d", (m %/% 60) %% 24, round(m %% 60)))

IN_BED <- c("in_bed", "asleep", "awake", "out_of_bed")

parse_night <- function(x) {
  tok <- trimws(strsplit(x, ",", fixed = TRUE)[[1]])
  val <- function(k) { i <- which(tok == k)[1]; if (is.na(i)) NA_character_ else tok[i + 1] }
  num <- function(k) suppressWarnings(as.numeric(val(k)))
  ti <- which(tok == "time")
  time <- tok[ti + 1]
  status <- tok[ti + 3]
  tot <- c(sl = num("sl"), wans = num("wans"), twt = num("twt"),
           tib = num("tib"), tst = num("tst"), se = num("se"))

  # first 96 epochs = 24 h from the grid's start (later epochs repeat the clock)
  s24 <- status[seq_len(min(96, length(status)))]
  nap <- 15 * sum(s24 == "nap")

  bed_idx  <- which(status %in% IN_BED)
  slp_idx  <- which(status == "asleep")
  final_wake <- ema <- n_wake <- out_night <- NA_real_
  split_night <- FALSE
  if (length(bed_idx)) {
    window <- min(bed_idx):max(bed_idx)
    split_night <- any(!status[window] %in% IN_BED)
    out_night <- 15 * sum(status[window] == "out_of_bed")
  }
  if (length(slp_idx)) {
    last <- max(slp_idx)
    final_wake <- to_min(time[last]) + 15
    after <- if (length(bed_idx) && max(bed_idx) > last) (last + 1):max(bed_idx) else integer(0)
    ema <- 15 * sum(status[after] %in% IN_BED)
    between <- status[min(slp_idx):last] %in% c("awake", "out_of_bed")
    r <- rle(between)
    n_wake <- sum(r$values)
  }
  data.frame(
    bed = val("start"), rise = val("end"),
    final_wake = fmt_clock(final_wake),
    sl = tot[["sl"]], wans = tot[["wans"]], twt = tot[["twt"]],
    tib = tot[["tib"]], tst = tot[["tst"]], se = tot[["se"]],
    ema = ema, n_wake = n_wake, out_night = out_night, nap = nap,
    split_night = split_night,
    grid_tst = 15 * sum(status == "asleep"),
    grid_tib = 15 * sum(status %in% IN_BED),
    stringsAsFactors = FALSE)
}

nights <- do.call(rbind, lapply(raw$Slaaptijd, parse_night))

# ---- checks: app totals must equal the grid ---------------------------------
chk <- c(
  "TST equals asleep epochs"   = all(nights$tst == nights$grid_tst, na.rm = TRUE),
  "TIB equals in-bed window"    = all(nights$tib == nights$grid_tib, na.rm = TRUE),
  "TIB = TST + SL + WANS"       = all(nights$tib == nights$tst + nights$sl + nights$wans, na.rm = TRUE),
  "EMA <= WANS"                 = all(is.na(nights$ema) | nights$ema <= nights$wans))
print(chk)
if (!all(chk)) warning("Some app totals disagree with the grid; see the 'controle' sheet.")

# ---- assemble one row per participant-night ---------------------------------
num <- function(x) suppressWarnings(as.numeric(gsub(",", ".", x)))
dag <- as.Date(num(raw$Dag), origin = "1899-12-30")

# Medication: free text. 0 = an explicit 'no', 1 = anything else written,
# NA = left empty. Check the 'medicatie_codering' sheet before using it.
med_txt <- trimws(raw$Medicatie)
# ("neem" alone is read as a typo of "neen")
med_no  <- tolower(med_txt) %in% c("nee", "neen", "neem", "non", "none", "false", "geen",
                                   "nvt", "n.v.t.", "no", "0", "-", "/", "niks", "niets",
                                   "geen medicatie")
med_bin <- ifelse(is.na(med_txt) | med_txt == "", NA_integer_, ifelse(med_no, 0L, 1L))

out <- data.frame(
  Owner_Id              = raw$`Owner Id`,
  Datum                 = dag,
  Naar_bed_gegaan       = nights$bed,
  Wakker_geworden       = nights$final_wake,
  Opgestaan             = nights$rise,
  Inslaaptijd_min       = nights$sl,
  Aantal_keer_wakker    = nights$n_wake,
  Hoelang_wakker_min    = nights$wans - nights$ema,
  Vroeg_wakker_min      = nights$ema,
  Buiten_bed_nacht_min  = nights$out_night,
  Totale_slaaptijd_min  = nights$tst,
  Tijd_in_bed_min       = nights$tib,
  Slaapefficientie      = nights$se,
  Dutjes_min            = nights$nap,
  Slaapkwaliteit        = num(raw$Slaapkwaliteit),
  Uitgerust             = num(raw$Uitgerust),
  Koffie                = num(raw$Koffie),
  Koffie_na_15u         = num(raw$`Koffie (na 15u)`),
  Alcohol               = num(raw$Alcohol),
  Tijd_afbouwen_min     = num(raw$`Tijd afbouwen`),
  Medicatie_ja          = med_bin,
  Medicatie             = med_txt,
  Extra_info            = raw$`Extra info`,
  Gesplitste_nacht      = nights$split_night,
  Ingevuld_op           = raw$`Created At`,
  Gewijzigd_op          = raw$`Updated At`,
  stringsAsFactors = FALSE)

# ---- clean ------------------------------------------------------------------
log <- character(0)
drop_no_id <- is.na(out$Owner_Id) | out$Owner_Id == ""
log <- c(log, sprintf("Rows without Owner Id removed: %d", sum(drop_no_id)))
removed <- out[drop_no_id, ]
removed$reden <- rep("geen Owner Id", nrow(removed))
out <- out[!drop_no_id, ]

# duplicate participant-days: keep the most recently updated entry
out <- out[order(out$Owner_Id, out$Datum, out$Gewijzigd_op, decreasing = c(FALSE, FALSE, TRUE),
                 method = "radix"), ]
dup <- duplicated(out[, c("Owner_Id", "Datum")])
log <- c(log, sprintf("Duplicate participant-days removed (latest edit kept): %d", sum(dup)))
if (any(dup)) {
  d2 <- out[dup, ]; d2$reden <- rep("dubbele dag (oudere versie)", nrow(d2))
  removed <- rbind(removed, d2)
}
out <- out[!dup, ]

# short participant codes in order of first entry
first <- tapply(out$Datum, out$Owner_Id, min)
codes <- names(sort(first))
key <- data.frame(Participant = sprintf("SVT%02d", seq_along(codes)), Owner_Id = codes,
                  stringsAsFactors = FALSE)
out$Participant <- key$Participant[match(out$Owner_Id, key$Owner_Id)]
out <- out[order(out$Participant, out$Datum), ]
out$Dag_nr <- as.integer(out$Datum - as.Date(first[out$Owner_Id]) + 1)
out <- out[, c("Participant", "Dag_nr", setdiff(names(out), c("Participant", "Dag_nr")))]
rownames(out) <- NULL

# flags (kept in the data; the Shiny app lists invalid nights itself)
flag <- character(nrow(out))
flag[out$Totale_slaaptijd_min == 0] <- "geen slaap geregistreerd"
flag[flag == "" & out$Tijd_in_bed_min > 16 * 60] <- "tijd in bed > 16 u"
flag[flag == "" & out$Gesplitste_nacht] <- "nacht onderbroken (uit bed en terug)"
out$Controle <- flag
log <- c(log, sprintf("Nights flagged '%s': %d", names(table(flag[flag != ""])),
                      as.integer(table(flag[flag != ""]))))
log <- c(log, sprintf("Result: %d nights from %d participants (median %d nights per participant)",
                      nrow(out), length(unique(out$Participant)),
                      as.integer(stats::median(table(out$Participant)))))

codebook <- data.frame(
  variabele = names(out),
  betekenis = c(
    "Korte deelnemerscode (volgorde van eerste invoer); koppeling met Owner_Id in blad 'sleutel'",
    "Dagnummer sinds de eerste invoer van deze deelnemer (gaten = ontbrekende nachten)",
    "Originele app-ID van de deelnemer",
    "Avond waarop de nacht begint ('Dag' in de export); de invoer gebeurt meestal de ochtend erna",
    "In bed gegaan (klok; 'start')",
    "Laatste ontwaken = einde van de laatste 15-min slaapperiode (klok)",
    "Uit bed (klok; 'end')",
    "Inslaaptijd (min; 'sl')",
    "Aantal wakkerperiodes tussen inslapen en laatste ontwaken (uit het 15-min rooster)",
    "WASO: wakker na inslapen tot het laatste ontwaken (min) = 'wans' - Vroeg_wakker_min",
    "EMA: wakker in bed tussen laatste ontwaken en opstaan (min)",
    "Tijd uit bed tijdens de nacht (min; status 'out_of_bed'), zit in WASO",
    "Totale slaaptijd (min; 'tst')",
    "Tijd in bed (min; 'tib' = TST + inslaaptijd + WASO + EMA)",
    "Slaapefficientie (%; 'se')",
    "Dutjes in de 24 u van het rooster (min; status 'nap')",
    "Slaapkwaliteit (0-10) over deze nacht",
    "Uitgerust (0-10) bij het opstaan na deze nacht",
    "Aantal koffie op 'Dag'",
    "Aantal koffie na 15 u op 'Dag'",
    "Aantal alcoholconsumpties op 'Dag'",
    "Tijd afbouwen voor het slapengaan (min)",
    "Medicatie genomen: 0 = expliciet nee, 1 = iets ingevuld, leeg = niet ingevuld (zie blad 'medicatie_codering')",
    "Medicatie (vrije tekst)",
    "Extra info (vrije tekst)",
    "TRUE als de bedperiode onderbroken werd door 'not_in_bed' of 'nap'",
    "Tijdstip van invoer in de app",
    "Tijdstip van laatste wijziging in de app",
    "Waarschuwing om na te kijken (leeg = geen)"),
  stringsAsFactors = FALSE)

med_tab <- as.data.frame(table(Medicatie = out$Medicatie, Medicatie_ja = out$Medicatie_ja,
                               useNA = "ifany"), stringsAsFactors = FALSE)
med_tab <- med_tab[med_tab$Freq > 0, ]
med_tab <- med_tab[order(-med_tab$Freq), ]

write_xlsx(list(
  nachten = out,
  codeboek = codebook,
  sleutel = key,
  medicatie_codering = med_tab,
  verwijderd = removed,
  log = data.frame(stap = c(names(chk), log),
                   resultaat = c(ifelse(chk, "OK", "AFWIJKING"), rep("", length(log))))
), out_file)

cat(log, sep = "\n")
cat("Written:", out_file, "\n")
