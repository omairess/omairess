# ==============================================================================
# server/02b_helpers_palette.R — ONE group palette for the whole app
#
# WHY THIS EXISTS
# ---------------
# Seven different palettes were in use: 'red'/'blue'/'green' in the fPCA tab,
# Set1 in FoSR, a Brewer-ish set in clustering, firebrick/steelblue/forestgreen
# in the harmonic fit plot, a pastel set in its group-comparison overlay, and
# FCK_DENSITY_COLORS in the polar density tab. The same age band was a different
# colour in every figure, which makes cross-reading them slow and error-prone.
#
# COLOUR FOLLOWS THE ENTITY, NEVER ITS RANK
# -----------------------------------------
# The subtler half of the problem: every one of those palettes was indexed by
# POSITION. Filter one group out, or sort the levels differently in one tab, and
# the survivors are repainted -- so "the blue group" means different things in
# two figures on the same screen. fck_group_colors() keys on the LEVEL NAME
# through a stable sorted level set, so a group keeps its colour no matter what
# is filtered, dropped or reordered.
#
# THE PALETTE ITSELF: "ALTERNATING LIGHT/DARK PAIRS"
# --------------------------------------------------
# The reporter supplied this palette, as four hue families each with a light and
# a medium member:
#
#     Light Blue  #8fd7d7   Med Blue   #00b0be
#     Light Pink  #ff8ca1   Med Pink   #f45f74
#     Light Green #bdd373   Med Green  #98c127
#     Light Orange#ffcd8e   Med Orange #ffb255
#
# ORDERING. Taken in the order listed, slots 1 and 2 would be the light and the
# medium BLUE -- two series that differ by OKLab dE 14.9 and read as one group
# split in two. A palette of this shape is built to be used one ROW at a time:
# the four hue families are the categorical dimension, light/dark is the
# alternation within a family. So the mediums take slots 1-4 and the lights take
# slots 5-8, which means slot k and slot k+4 are the two members of one family
# and the common 2-, 3- and 4-group figures get the four maximally separated
# hues.
#
# WHAT THE VALIDATOR SAYS ABOUT IT (light mode, --pairs all, the right mode
# here: a scatter or an overlay compares every pair, not just neighbours)
#
#   4 slots, mediums (the reporter's data has four age bands):
#     [FAIL] Lightness band    #ffb255 L 0.820, above the 0.77 ceiling
#     [PASS] Chroma floor
#     [FAIL] CVD separation    #ffb255 <-> #98c127  dE 3.2 (protan)
#     [PASS] Normal-vision     worst dE 16.5
#     [WARN] Contrast          3 of 4 below 3:1 vs the surface
#
#   8 slots: adds a chroma failure (#8fd7d7 C 0.072 and #ffcd8e C 0.098 read as
#   gray), and every light/medium pair within a family is below the 15 dE
#   normal-vision floor -- 7.2 for the two oranges. That is inherent to the
#   design: the pairs are MEANT to read as related. It is not a defect of the
#   palette, it is what the palette is for. Past four groups this app therefore
#   carries identity in LINE DASH as well as hue, and the figures say so.
#
# THE ONE THING THAT MATTERS FOR THIS DATASET is the medium green against the
# medium orange: dE 3.2 under protanopia means two of the four age bands are the
# same colour for roughly 1 in 12 male readers. FCK_PALETTE_VARIANT switches
# between:
#
#   "supplied" (default) -- the eight hexes exactly as given.
#   "cvd_safe"           -- the same four hue families with the orange deepened
#                           (#ffb255 -> #d38400, same OKLCH hue angle 68.6, L
#                           0.82 -> 0.68) and the green nudged (#98c127 ->
#                           #9ec055). That clears every hard check: lightness
#                           PASS, chroma PASS, normal-vision PASS (worst 15.2),
#                           CVD at the 8.0 target boundary (WARN, relieved by
#                           the legend and the tables every one of these figures
#                           carries).
#
# Change the one constant below to switch; nothing else in the app needs to
# know. The app is light-mode only (shinydashboard on a light surface), so the
# light validation is the applicable one.
# ==============================================================================

# "supplied" = the reporter's hexes verbatim; "cvd_safe" = green/orange
# re-stepped so all-pairs CVD separation clears. See the header.
FCK_PALETTE_VARIANT <- "supplied"

# The four hue families, medium member -- slots 1-4.
.FCK_MEDIUMS <- list(
  supplied = c("#00b0be", "#f45f74", "#98c127", "#ffb255"),
  cvd_safe = c("#00b0be", "#f45f74", "#9ec055", "#d38400")
)

# The same four hue families, light member -- slots 5-8. Only ever used
# alongside a dash pattern; see FCK_GROUP_DASH.
.FCK_LIGHTS <- c("#8fd7d7", "#ff8ca1", "#bdd373", "#ffcd8e")

FCK_GROUP_COLORS <- .FCK_MEDIUMS[[
  if (FCK_PALETTE_VARIANT %in% names(.FCK_MEDIUMS)) FCK_PALETTE_VARIANT else "supplied"
]]

# Hues past the validated four. They are NOT all-pairs separable from the first
# four -- each is the light partner of one of them -- and are only ever used
# together with a dash pattern.
FCK_GROUP_COLORS_EXTRA <- .FCK_LIGHTS

# Beyond four slots hue alone is not enough; dash carries identity with it.
# Eight distinct patterns, so slot k and its light partner at slot k+4 never
# share one.
FCK_GROUP_DASH <- c("solid", "dash", "dot", "dashdot",
                    "longdash", "longdashdot", "dash", "dot")

# The full ramp a caller may draw from, validated prefix first.
fck_group_ramp <- function(n) {
  full <- c(FCK_GROUP_COLORS, FCK_GROUP_COLORS_EXTRA)
  if (n <= length(full)) full[seq_len(n)] else rep_len(full, n)
}

# Colours for the harmonic COMPONENTS (H1, H2, ...) and for principal
# components -- a different entity class from groups, but drawn from the same
# palette so one screen never mixes two colour vocabularies. These replaced
# c("green","orange","purple","brown","pink","cyan","magenta","olive"), which
# were unvalidated CSS names.
fck_component_colors <- function(n) fck_group_ramp(n)

# ---- the one function every plot should call --------------------------------
#
# Give it the group levels (in any order, with or without NAs) and it returns a
# NAMED colour vector. Look up by name, never by position:
#
#     cols <- fck_group_colors(levels(g))
#     cols[[as.character(this_group)]]
#
# `all_levels` pins the assignment to the full level set even when the current
# figure shows a subset, so filtering one group out never repaints the others.
fck_group_colors <- function(levels, all_levels = NULL) {
  lv <- as.character(levels)
  lv <- lv[!is.na(lv) & nzchar(lv)]
  ref <- if (!is.null(all_levels)) {
    r <- as.character(all_levels); r[!is.na(r) & nzchar(r)]
  } else lv
  ref <- sort(unique(ref))   # the same order fck_group_reference() produces
  pal <- fck_group_ramp(length(ref))
  names(pal) <- ref
  out <- pal[unique(lv)]
  names(out) <- unique(lv)
  out[is.na(out)] <- "#6f6f6f"   # a level not in the reference set reads gray
  out
}

# The reference order both the colour and the dash lookups index into. Sorting
# is what makes two tabs that received the levels in different orders agree.
fck_group_reference <- function(levels, all_levels = NULL) {
  ref <- as.character(if (is.null(all_levels)) levels else all_levels)
  ref <- ref[!is.na(ref) & nzchar(ref)]
  sort(unique(ref))
}

# The matching dash pattern, keyed exactly as the colours are.
#
# This originally took its reference order from names(fck_group_colors(...)),
# which are in INPUT order rather than the sorted order the colours are assigned
# from -- so a group's dash changed when the levels arrived differently even
# though its colour did not. Both now index the same reference.
fck_group_dashes <- function(levels, all_levels = NULL) {
  ref <- fck_group_reference(levels, all_levels)
  lv <- unique(as.character(levels))
  lv <- lv[!is.na(lv) & nzchar(lv)]
  idx <- match(lv, ref)
  d <- FCK_GROUP_DASH[((idx - 1) %% length(FCK_GROUP_DASH)) + 1]
  d[is.na(idx)] <- "solid"
  names(d) <- lv
  d
}

# Same colour with an alpha suffix, for confidence bands and fills. plotly takes
# 8-digit hex, so a band and its line cannot drift apart.
fck_group_fill <- function(hex, alpha = 0.25) {
  a <- sprintf("%02X", as.integer(round(pmax(0, pmin(1, alpha)) * 255)))
  paste0(hex, a)
}

# rgba(), for the places that build a colour string by hand.
fck_group_rgba <- function(hex, alpha = 0.25) {
  m <- grDevices::col2rgb(hex)
  sprintf("rgba(%d,%d,%d,%.3f)", m[1, ], m[2, ], m[3, ], alpha)
}

# Does this figure need dash as well as hue to stay readable? True past the four
# hue families, where the fifth slot is the light partner of the first.
fck_needs_dash <- function(n) n > length(FCK_GROUP_COLORS)

# ---- roles that are NOT groups ----------------------------------------------
#
# A figure that draws one summary line over one data series still needs colour,
# and the harmonic tab used 'red', 'firebrick' and 'steelblue' for it. Red in
# particular is a problem now: the palette's slot 2 is a pink (#f45f74), so a
# red population-mean line reads as a fifth group. These three roles are
# deliberately OUTSIDE the categorical palette.

# The summary line -- a population mean, a pooled fit. Dark neutral: it is the
# emphasis, not a category, and it must never be confused for one.
FCK_EMPHASIS <- "#33322e"

# Reference geometry that carries no data: identity lines, y = 0, a confidence
# ellipse over the whole sample.
FCK_NEUTRAL <- "#7a7873"

# The single data series in a figure that has no grouping. Slot 1, so a plot
# that later gains a grouping variable does not change colour vocabulary.
FCK_SERIES1 <- unname(FCK_GROUP_COLORS[1])
FCK_SERIES2 <- unname(FCK_GROUP_COLORS[2])

# Status TEXT on a figure ("the fit did not converge"), not a data series.
# Deliberately a dark red rather than plotly's 'red', which at the palette's
# pink is close enough to read as a category.
FCK_ALERT <- "#a3261f"
