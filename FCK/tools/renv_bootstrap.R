# ==============================================================================
# tools/renv_bootstrap.R — record the environment this app runs in
# ==============================================================================
# AUDIT (P2.2). app.R used to call install.packages() while starting; that is now
# removed. It made the analysis depend on whatever CRAN held on the day, and it
# could not work on the deployed or offline machines where it mattered most.
#
# The replacement is a lockfile. It is deliberately NOT committed by the audit:
# a lockfile records the library an analysis actually ran against, and the audit
# container is not that machine (fda.usc, refund, reticulate and shinyWidgets
# are not installed there, so a lockfile written from it would be a fiction).
# Run this on the machine you analyse on, and commit what it produces.
#
#   Rscript tools/renv_bootstrap.R          # from the FCK directory
#
# Then, to restore that exact library later or elsewhere:
#
#   renv::restore()
# ==============================================================================

if (!requireNamespace("renv", quietly = TRUE)) {
  stop("renv is not installed. Run:  install.packages(\"renv\")\n",
       "This script deliberately does not install it for you: the reason the\n",
       "startup installer was removed is that a tool should not silently change\n",
       "the library your results depend on.", call. = FALSE)
}

app_root <- getwd()
if (!dir.exists(file.path(app_root, "server")))
  stop("Run this from the FCK directory (the one containing app.R and server/).",
       call. = FALSE)
message("Project root: ", app_root)

# renv finds dependencies by scanning the source, so every library() and pkg::
# call in the app is recorded -- including the optional ones, which is the point.
renv::init(project = app_root, bare = FALSE, restart = FALSE)
renv::snapshot(project = app_root, prompt = FALSE)

message("\nWritten: ", file.path(app_root, "renv.lock"))
message("Commit renv.lock alongside your analysis.\n")

cat("\nFor the record, this machine is running:\n")
cat(" ", R.version.string, "\n")
for (p in c("shiny", "fda", "mgcv", "refund", "fda.usc", "minpack.lm",
            "plotly", "ggplot2", "cluster", "readxl")) {
  v <- tryCatch(as.character(utils::packageVersion(p)), error = function(e) NA_character_)
  cat(sprintf("  %-12s %s\n", p, if (is.na(v)) "not installed" else v))
}
