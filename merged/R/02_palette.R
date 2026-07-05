# ============================================================================
# 02_palette.R — pastel defaults (house rule 9)
#
# Single source of default colours for all three tabs. Sourced from
# RColorBrewer Pastel1/Pastel2 as the rule requires — not the three ad hoc
# palettes the source apps shipped (rainbow(), "Blues" sequential, hardcoded
# hex arrays), all of which are dropped.
# ============================================================================

# DECISION: positive edge = Pastel1 green, negative edge = Pastel1 red-pink,
#   matching the qgraph green/red sign convention but in pastel; node fill =
#   Pastel1 blue. Dagger's arcs are UNSIGNED, so its tab maps "positive edge
#   colour" onto arc colour and ignores the negative picker (documented in
#   the appearance module) rather than pretending bnlearn arcs have sign.
house_pastel <- function() {
  p1 <- RColorBrewer::brewer.pal(9, "Pastel1")
  list(
    node_fill   = p1[2],   # "#B3CDE3" light blue
    node_border = "#7F7F7F",
    pos_edge    = p1[3],   # "#CCEBC5" light green
    neg_edge    = p1[1],   # "#FBB4AE" light red-pink
    groups      = p1       # categorical fills (communities, factors)
  )
}

# Categorical fills for n groups; interpolate past Brewer's 9-colour maximum.
house_group_colors <- function(n) {
  p1 <- RColorBrewer::brewer.pal(9, "Pastel1")
  if (n <= 9) p1[seq_len(n)]
  else grDevices::colorRampPalette(p1)(n)
}
