# ============================================================================
# 02_palette.R — default colours
#
# Single source of default colours for all tabs. The user-chosen defaults:
#   nodes #FFFFFF, node border #878787, positive edges #2C4778 (deep blue),
#   negative edges #E03920 (red), predictability ring #FFAE00 (amber).
#
# For CATEGORICAL fills (communities, latent factors) we use a custom
# 12-colour qualitative palette; past that count colours are interpolated.
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

HOUSE_QUAL_PALETTE <- c(
  "#4F97BC", # muted blue
  "#FDB24A", # soft orange
  "#45B8BC", # turquoise
  "#DFAF47", # ochre
  "#E15F39", # soft vermillion
  "#F4C21A", # warm yellow
  "#EFD96A", # pale yellow
  "#C47A35", # burnt orange
  "#B8CDBF", # sage
  "#C8C481", # khaki
  "#697A20", # olive green
  "#7A4A18"  # brown
)

# Categorical fills for n groups; interpolate past the curated set for large n.
house_group_colors <- function(n) {
  if (n <= length(HOUSE_QUAL_PALETTE)) HOUSE_QUAL_PALETTE[seq_len(n)]
  else grDevices::colorRampPalette(HOUSE_QUAL_PALETTE)(n)
}
