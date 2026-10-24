# ============================================================================
# 00_packages.R — package guard + version ledger (house rule 0)
#
# Installs only what is MISSING; never updates an installed package (a silent
# brms/Stan update can change posterior summaries). Records exact versions for
# the script header and the export bundle.
# ============================================================================

APP_PACKAGES <- c(
  # UI
  "shiny", "shinyWidgets", "colourpicker", "DT", "svglite",
  # import
  "readr", "readxl", "haven", "jsonlite", "tidyr", "tidyselect", "stringi",
  "data.table", "writexl",
  # analysis + figures
  "multilevelcoda", "brms", "compositions", "ggplot2", "RColorBrewer"
)

# Packages the exported script needs (no Shiny).
SCRIPT_PACKAGES <- c("readr", "readxl", "haven", "jsonlite", "tidyr",
                     "tidyselect", "stringi", "data.table", "ggplot2",
                     "RColorBrewer", "brms", "multilevelcoda")

# Install missing packages only. Updating is never done here.
ensure_packages <- function(pkgs = APP_PACKAGES) {
  missing <- setdiff(pkgs, rownames(utils::installed.packages()))
  if (length(missing)) {
    utils::install.packages(missing, repos = "https://cloud.r-project.org")
  }
  # Loaded in this order so multilevelcoda's methods win over brms' generics;
  # multilevelcoda masks base::sub(), which the helpers never call unqualified.
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

# Is a C++ toolchain available for Stan? brms compiles every model.
stan_toolchain_ok <- function() {
  ok <- tryCatch(pkgbuild::has_build_tools(debug = FALSE),
                 error = function(e) NA)
  isTRUE(ok) || is.na(ok)
}
