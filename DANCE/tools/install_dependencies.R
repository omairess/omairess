#!/usr/bin/env Rscript
# ==============================================================================
# tools/install_dependencies.R — install everything DANCE can use, in one go
# ==============================================================================
# The app itself asks before installing, and only in an interactive session.
# This is the non-interactive equivalent: for a fresh machine, a container
# build, or anyone who would rather run one command than answer a prompt.
#
#   Rscript tools/install_dependencies.R            # required + optional
#   Rscript tools/install_dependencies.R --required # the app's hard minimum
#
# It reads the two package vectors OUT OF app.R rather than repeating them, so
# this file cannot drift from what the app actually checks for.
# ==============================================================================

args <- commandArgs(trailingOnly = TRUE)
only_required <- "--required" %in% args
root <- if (file.exists("app.R")) "." else if (file.exists("DANCE/app.R")) "DANCE" else
  stop("Run this from the DANCE directory (or its parent).", call. = FALSE)

env <- new.env()
for (ex in parse(file.path(root, "app.R"))) {
  if (is.call(ex) && identical(as.character(ex[[1]]), "<-") && is.name(ex[[2]]) &&
      as.character(ex[[2]]) %in% c("required_packages", "optional_packages"))
    eval(ex, env)
}
want <- c(env$required_packages, if (!only_required) names(env$optional_packages))
want <- unique(want)
have <- vapply(want, requireNamespace, logical(1), quietly = TRUE)

cat(sprintf("DANCE dependencies: %d wanted, %d already installed\n", length(want), sum(have)))
if (!any(!have)) { cat("Nothing to do.\n"); quit(status = 0) }

cat("Installing:", paste(want[!have], collapse = ", "), "\n\n")
install.packages(want[!have])

still <- want[!vapply(want, requireNamespace, logical(1), quietly = TRUE)]
cat("\n")
for (p in setdiff(want, still))
  cat(sprintf("  %-14s %s\n", p, as.character(packageVersion(p))))
if (length(still)) {
  cat("\nSTILL MISSING:", paste(still, collapse = ", "), "\n")
  # Optional ones are survivable and the app says what each costs; a missing
  # required one means the app will not start, so the exit status differs.
  req_missing <- intersect(still, env$required_packages)
  if (length(req_missing)) {
    cat("These are REQUIRED -- the app will not start without them.\n")
    quit(status = 1)
  }
  cat("These are optional; the app will start and report what they disable.\n")
}
