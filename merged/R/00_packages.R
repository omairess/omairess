# ============================================================================
# 00_packages.R — package guard + version ledger (house rule 0)
#
# Installs only what is MISSING; never updates an installed package (a silent
# qgraph/bootnet/glasso update can change a reported network). Records exact
# versions for the script header and the export bundle.
# ============================================================================

APP_PACKAGES <- c(
  # shared UI
  "shiny", "shinyWidgets", "colourpicker", "DT", "svglite",
  # shared data handling
  "readr", "readxl", "haven", "tidyr", "tools", "RColorBrewer",
  # bootnet tab (estimators: EBICglasso/ggmModSelect/pcor/cor/huge/
  #   IsingFit/mgm/relimp; EGAnet for the EGA community layout)
  "bootnet", "qgraph", "networktools", "huge", "IsingFit", "mgm",
  "relaimpo", "EGAnet",
  # DAG tab
  "bnlearn", "igraph",
  # psychonetrics / NCT tab
  "psychonetrics", "NetworkComparisonTest", "semPlot"
)

# Install missing packages only. Updating is opt-in elsewhere, never here.
ensure_packages <- function(pkgs = APP_PACKAGES) {
  missing <- setdiff(pkgs, rownames(installed.packages()))
  if (length(missing)) {
    install.packages(missing, repos = "https://cloud.r-project.org")
  }
  for (p in pkgs) suppressPackageStartupMessages(
    library(p, character.only = TRUE)
  )
  invisible(pkgs)
}

# Exact versions, for the script header (rule 6) and export bundle (rule 7).
pkg_versions <- function(pkgs = APP_PACKAGES) {
  vapply(pkgs, function(p) {
    tryCatch(as.character(utils::packageVersion(p)),
             error = function(e) "NOT INSTALLED")
  }, character(1))
}

# Record the version ledger as the first pipeline step, so every exported
# script and log opens with the environment it was produced in.
record_setup_step <- function(rec, pkgs = APP_PACKAGES) {
  vers <- pkg_versions(pkgs)
  rec_upsert(
    rec, id = "shared_setup", phase = "setup",
    description = sprintf(
      "Loaded %d packages (versions recorded in the script header); no automatic package updates.",
      length(pkgs)),
    code = paste(
      "# Versions this analysis was produced with are in the header above.",
      "# Packages were installed only if missing — never auto-updated.",
      sep = "\n")
  )
  vers
}
