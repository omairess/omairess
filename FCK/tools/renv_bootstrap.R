# ==============================================================================
# tools/renv_bootstrap.R — record the environment this app runs in
# ==============================================================================
# AUDIT (P2.2, extended at P3.6). app.R used to call install.packages() while
# starting; that is now removed. It made the analysis depend on whatever CRAN
# held on the day, and it could not work on the deployed or offline machines
# where it mattered most. The replacement is a lockfile.
#
# WHY renv.lock IS NOT COMMITTED IN THIS REPOSITORY
#
# A reviewer asked for one to be committed. It is not, and this is the reason,
# stated plainly so the next person does not have to re-derive it.
#
# A lockfile is a RECORD OF A LIBRARY THAT EXISTS. renv::snapshot() writes the
# version and hash of each package as installed on the machine it runs on. The
# container this app was audited in does not have fda.usc, shinyWidgets,
# reticulate or tidyr installed. A lockfile written there would either omit
# them -- in which case renv::restore() silently does not install them, and the
# clustering and cosinor tabs fail at first use with a missing-package error on
# a machine whose owner was told the environment was pinned -- or carry
# invented versions and hashes, which is worse: restore() would fail on the
# hash, or install something nobody verified.
#
# A lockfile that does not describe a working library is not a weaker
# guarantee than none. It is a false one, and this app is being corrected
# precisely because it made claims of that shape (a "REML" comment over an
# unpenalised fit, an "rmfanova" branch that never called rmfanova, an
# "exported analysis code" that could not run). Committing a fabricated
# lockfile would be another.
#
# So: run this on the machine you analyse on, and commit what it produces.
# The check below refuses to write a lockfile that is missing a REQUIRED
# package, which is exactly what stops the audit container from producing one.
#
#   Rscript tools/renv_bootstrap.R          # from the FCK directory
#
# Then, to restore that exact library later or elsewhere:
#
#   renv::restore()
# ==============================================================================

app_root <- getwd()
if (!dir.exists(file.path(app_root, "server")) ||
    !file.exists(file.path(app_root, "app.R")))
  stop("Run this from the FCK directory (the one containing app.R and server/).",
       call. = FALSE)

# ---- The package lists come from app.R, not from a copy kept here -----------
# P3.6: a second hand-maintained list is a second thing that goes stale. app.R
# is the authority on what this app needs; read its declarations directly.
app_env <- new.env()
app_exprs <- parse(file.path(app_root, "app.R"), encoding = "UTF-8")
for (ex in app_exprs) {
  if (is.call(ex) && length(ex) == 3 && as.character(ex[[1]]) %in% c("<-", "=") &&
      is.name(ex[[2]]) &&
      as.character(ex[[2]]) %in% c("required_packages", "optional_packages"))
    eval(ex, envir = app_env)
}
required <- get0("required_packages", envir = app_env)
optional <- get0("optional_packages", envir = app_env)
if (is.null(required))
  stop("Could not read required_packages from app.R. Has that block moved?",
       call. = FALSE)

have <- function(p) requireNamespace(p, quietly = TRUE)
ver  <- function(p) tryCatch(as.character(utils::packageVersion(p)),
                             error = function(e) NA_character_)

message("Project root: ", app_root)
message("\nRequired packages (", length(required), "):")
for (p in required)
  message(sprintf("  %-14s %s", p, if (have(p)) ver(p) else "MISSING"))

message("\nOptional packages (", length(optional), "):")
for (i in seq_along(optional))
  message(sprintf("  %-14s %-10s  %s", names(optional)[i],
                  if (have(names(optional)[i])) ver(names(optional)[i]) else "not installed",
                  optional[[i]]))

missing_req <- required[!vapply(required, have, logical(1))]
if (length(missing_req))
  stop("\nNot writing a lockfile: these REQUIRED packages are not installed here.\n",
       "  ", paste(missing_req, collapse = ", "), "\n\n",
       "A lockfile records a library that exists. One written from a machine\n",
       "missing these would not describe a working environment, and\n",
       "renv::restore() would reproduce the gap rather than fill it.\n",
       "Install them first:\n",
       "  install.packages(c(", paste(sprintf('"%s"', missing_req), collapse = ", "), "))",
       call. = FALSE)

missing_opt <- names(optional)[!vapply(names(optional), have, logical(1))]
if (length(missing_opt)) {
  message("\nNOTE: these optional packages are not installed, so they will NOT")
  message("appear in the lockfile, and renv::restore() will not bring them back:")
  for (p in missing_opt)
    message(sprintf("  %-14s (%s)", p, optional[[p]]))
  message("Install any feature you intend to use before snapshotting.")
}

if (!requireNamespace("renv", quietly = TRUE)) {
  stop("\nrenv is not installed. Run:  install.packages(\"renv\")\n",
       "This script deliberately does not install it for you: the reason the\n",
       "startup installer was removed is that a tool should not silently change\n",
       "the library your results depend on.", call. = FALSE)
}

# renv finds dependencies by scanning the source, so every library() and pkg::
# call in the app is recorded -- including the optional ones, which is the point.
renv::init(project = app_root, bare = FALSE, restart = FALSE)
renv::snapshot(project = app_root, prompt = FALSE)

message("\nWritten: ", file.path(app_root, "renv.lock"))
message("Commit renv.lock alongside your analysis.\n")

cat("\nFor the record, this machine is running:\n")
cat(" ", R.version.string, "\n")
for (p in c(required, names(optional)))
  cat(sprintf("  %-14s %s\n", p, if (is.na(ver(p))) "not installed" else ver(p)))
