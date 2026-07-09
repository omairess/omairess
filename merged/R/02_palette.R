# ============================================================================
# 02_palette.R — default colours
#
# Single source of default colours for all tabs. The user-chosen defaults:
#   nodes #FFFFFF, node border #878787, positive edges #2C4778 (deep blue),
#   negative edges #E03920 (red), predictability ring #FFAE00 (amber).
#
# For CATEGORICAL fills (communities, latent factors) we drop the pastel set
# in favour of a curated qualitative palette following the guidance in
# datawrapper.de/blog/beautifulcolors — colours that are perceptually
# DISTINCT, not over-saturated, and reasonably colour-vision-deficiency
# friendly. This is the Okabe-Ito 8-colour set (a widely used CVD-safe
# qualitative palette) plus a few extra distinct hues; > that count is
# interpolated. Okabe-Ito is the standard example of the "distinguishable,
# not garish" principle the article advocates.
# ============================================================================

house_pastel <- function() {
  list(
    node_fill   = "#FFFFFF",
    node_border = "#878787",
    pos_edge    = "#2C4778",   # deep blue
    neg_edge    = "#E03920",   # red
    ring        = "#FFAE00",   # amber (predictability rings)
    groups      = HOUSE_QUAL_PALETTE
  )
}

# Curated qualitative palette (Okabe-Ito + 4 extra distinct hues). Ordered so
# the first few are maximally distinct — good for the common 2-4 community /
# factor case.
HOUSE_QUAL_PALETTE <- c(
  "#2C4778", # deep blue   (matches the positive-edge colour)
  "#E69F00", # orange
  "#009E73", # green
  "#CC79A7", # pink
  "#56B4E9", # sky blue
  "#D55E00", # vermillion
  "#F0E442", # yellow
  "#7B3294", # purple
  "#117733", # dark green
  "#882255", # wine
  "#44AA99", # teal
  "#999933"  # olive
)

# Categorical fills for n groups; interpolate past the curated set for large n.
house_group_colors <- function(n) {
  if (n <= length(HOUSE_QUAL_PALETTE)) HOUSE_QUAL_PALETTE[seq_len(n)]
  else grDevices::colorRampPalette(HOUSE_QUAL_PALETTE)(n)
}
