# ==============================================================================
# ui/00_theme.R — shared styling
#
# Union of the two source apps' <style> blocks (WaPaa1_3.R lines 63-77 and
# CIRCAREG.R lines 56-61); the class names are referenced from ported tabs of
# both families, so both sets are kept.
# ==============================================================================

ui_theme_css <- tags$head(
  tags$style(HTML("
    .content-wrapper, .right-side {
      background-color: #f4f4f4;
    }
    /* from WaPaa: functional ANOVA / post-hoc tabs */
    .significance-legend {
      background-color: #f0f0f0;
      padding: 10px;
      border-radius: 5px;
      margin: 10px 0;
    }
    .pairwise-summary-box {
      background-color: #e8f4fd;
      padding: 15px;
      border-radius: 8px;
      margin: 10px 0;
      border-left: 4px solid #2196F3;
    }
    /* from CIRCAREG: regression tabs */
    .reg-control-panel {
      background-color: #fff;
      padding: 15px;
      border: 1px solid #ddd;
      border-radius: 5px;
      margin-bottom: 15px;
    }
    .box-header { font-weight: bold; }
  "))
)
