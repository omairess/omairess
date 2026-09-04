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
# THE PALETTE IS VALIDATED, NOT CHOSEN BY EYE
# -------------------------------------------
# Checked with the dataviz validator in --pairs all mode, which is the right
# mode here: a scatter or an overlay compares every pair, not just neighbours.
# The previous palette FAILED that check at four groups (amber vs orange, normal
# -vision dE 13.7 -- below the 15 floor), which matters directly: the reporter's
# data has four age bands. Slot 4 was re-stepped to magenta and slot 5 to a
# saturated brown.
#
#   light mode, --pairs all, 5 slots: 0 failures
#     worst CVD separation        dE 9.2  (green vs orange, deutan)
#     worst normal-vision         dE 19.1
#     contrast WARN on the green  -> relief is the legend and the tables, which
#                                    every one of these figures carries
#
# FIVE is the honest maximum. A sixth all-pairs-distinct hue does not exist
# inside the lightness band -- every candidate tried failed the normal-vision
# floor. Past five, identity is carried by LINE DASH as well as hue, and the
# figures say so. A ninth series is not a generated hue; group the tail.
#
# The app is light-mode only (shinydashboard on a light surface), so the light
# validation is the applicable one; the orange sits outside the dark-mode
# lightness band and would need re-stepping if a dark theme is ever added.
# ==============================================================================

FCK_GROUP_COLORS <- c(
  "#2a78d6",  # blue
  "#eb6834",  # orange
  "#1baf7a",  # green
  "#b5309b",  # magenta
  "#8a5a12"   # brown
)

# Beyond five slots hue alone is not enough; dash carries identity with it.
FCK_GROUP_DASH <- c("solid", "dash", "dot", "dashdot", "longdash",
                    "solid", "dash", "dot")

# Hues past the validated five. They are NOT all-pairs separable from the first
# five and are only ever used together with a dash pattern.
FCK_GROUP_COLORS_EXTRA <- c("#4a3aa7", "#00707d", "#c2185b")

# The full ramp a caller may draw from, validated prefix first.
fck_group_ramp <- function(n) {
  full <- c(FCK_GROUP_COLORS, FCK_GROUP_COLORS_EXTRA)
  if (n <= length(full)) full[seq_len(n)] else rep_len(full, n)
}

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

# Does this figure need dash as well as hue to stay readable?
fck_needs_dash <- function(n) n > length(FCK_GROUP_COLORS)
