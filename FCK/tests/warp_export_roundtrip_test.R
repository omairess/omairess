# ==============================================================================
# tests/warp_export_roundtrip_test.R — does the EXPORTED script reproduce the
# registration the app actually ran?
#
# WHY THIS EXISTS. Every other guarantee in this audit is about the live app.
# This one is about the promise printed at the top of the exported script: that
# running it reproduces the analysis. For time warping that promise was false
# for all three methods, and had been since before P0.8, because
# server/90_export.R wrote its own copy of the algorithms with add("...") and
# the copy was never updated when the estimators were fixed. Measured then:
#
#     true delta    live s     exported s
#      0.05        -0.0505       0.0050      <- a twelfth of it, sign reversed
#      0.10        -0.1010       0.0090
#
# and on curves needing NO parametric warping the exported quadratic family
# deformed the time axis by 0.4999 of the domain where the live kernel gave
# 0.0003. tests/codegen_test.R could not catch it: an exported script full of
# the WRONG algorithm parses perfectly.
#
# The fix was structural (P11.1): the kernels moved to server/05_helpers_warp.R
# and the generator emits them with deparse(). This file is what stops that
# regressing. It does not compare the app against a re-implementation — it runs
# the SHIPPED kernel, and separately runs the generator's emitted text, and
# requires the two to agree to 1e-10.
#
# Run with:   Rscript tests/warp_export_roundtrip_test.R   (from the FCK directory)
# ==============================================================================
.libPaths(c("~/Rlib", .libPaths()))
ok <- 0L; bad <- 0L
chk <- function(cond, good, bad_msg) {
  if (isTRUE(cond)) { cat("ok   ", good, "\n"); ok <<- ok + 1L }
  else { cat("FAIL:", bad_msg, "\n"); bad <<- bad + 1L }
}

need <- c("fda", "codetools")
miss <- need[!vapply(need, requireNamespace, logical(1), quietly = TRUE)]
if (length(miss)) { cat("SKIP: missing", paste(miss, collapse = ", "), "\n"); quit(status = 0) }
suppressMessages(library(fda))

kernels <- c("linear_shift_alignment", "parametric_alignment",
             "landmark_alignment_simple", "fck_landmark_warp")

# ---- 1. the shipped kernel ---------------------------------------------------
kern <- new.env(parent = globalenv())
source("server/05_helpers_warp.R", local = kern)
for (f in kernels)
  chk(exists(f, envir = kern, inherits = FALSE),
      paste(f, "is a top-level kernel"),
      paste(f, "is NOT reachable at top level"))

# P11.2: the caller must be able to SEE its helper. This is the defect that made
# every landmark registration fall back to a linear shift, silently.
chk(identical(environment(kern$landmark_alignment_simple),
              environment(kern$fck_landmark_warp)),
    "landmark_alignment_simple and fck_landmark_warp share an environment",
    "fck_landmark_warp is not in its caller's scope")

# ---- 2. the kernels must be PURE, or they cannot be exported -----------------
for (f in kernels) {
  g <- codetools::findGlobals(get(f, envir = kern), merge = TRUE)
  reactive <- intersect(g, c("values", "input", "session", "showNotification",
                             "reactive", "observeEvent", "req"))
  chk(length(reactive) == 0,
      paste(f, "is pure"),
      paste0(f, " touches the reactive world (", paste(reactive, collapse = ", "),
             "), so deparse() cannot export it"))
}

# ---- 3. data with KNOWN periodic displacements ------------------------------
tp <- seq(0, 1, length.out = 100)
prof <- function(t) sin(2*pi*t) + 0.4*sin(4*pi*t + 0.8)
deltas <- c(0, 0.05, 0.10, -0.05, -0.10)
Y  <- sapply(deltas, function(d) prof(tp - d))
fd <- smooth.basis(tp, Y, fdPar(create.fourier.basis(c(0,1), nbasis = 15), 2, 1e-8))$fd

# ---- 4. run the generator's emitted text, and compare -----------------------
# The generator resolves kernels out of environment(generate_analysis_code); the
# kernel env stands in for the server environment here. What is under test is
# that emit_kernel() writes a definition which RUNS and AGREES — not that some
# hand-written copy happens to look similar.
emitted <- character(0)
add <- function(...) emitted <<- c(emitted, paste0(...))
emitted_fns <- character(0)
kernel_env  <- kern
fck_is_app_fn <- function(nm) {
  o <- tryCatch(get0(nm, envir = kernel_env, inherits = TRUE), error = function(e) NULL)
  is.function(o) && !is.null(environment(o)) && identical(environment(o), kernel_env)
}
fck_fn_order <- function(nm, seen = character(0), stack = character(0)) {
  if (nm %in% seen || nm %in% stack) return(seen)
  stack <- c(stack, nm)
  obj <- tryCatch(get(nm, envir = kernel_env), error = function(e) NULL)
  if (!is.function(obj)) return(seen)
  g <- tryCatch(codetools::findGlobals(obj, merge = TRUE), error = function(e) character(0))
  for (d in Filter(fck_is_app_fn, g)) seen <- fck_fn_order(d, seen, stack)
  c(seen, nm)
}
emit_kernel <- function(nm) {
  for (f in fck_fn_order(nm)) {
    if (f %in% emitted_fns) next
    obj <- tryCatch(get(f, envir = kernel_env), error = function(e) NULL)
    if (!is.function(obj)) next
    emitted_fns <<- c(emitted_fns, f)
    lhs <- if (identical(f, make.names(f))) f else paste0("`", f, "`")
    add(paste(lhs, "<-", paste(deparse(obj), collapse = "\n")))
    add("")
  }
}
run_emitted <- function(fn_name, call_txt) {
  emitted <<- character(0); emitted_fns <<- character(0)
  emit_kernel(fn_name)
  e <- new.env(parent = globalenv())
  assign("fd_obj", fd, envir = e); assign("time_points", tp, envir = e)
  eval(parse(text = paste(c(emitted, call_txt), collapse = "\n")), envir = e)
  get("warp_result", envir = e)
}

# -- linear shift, periodic (the case P10.1 corrected) ------------------------
live <- kern$linear_shift_alignment(fd, periodic = TRUE, reference = "first",
                                    time_points = tp)
expd <- run_emitted("linear_shift_alignment",
  "warp_result <- linear_shift_alignment(fd_obj, periodic = TRUE, reference = 'first', time_points = time_points)")
d1 <- max(abs(live$shifts - expd$shifts))
chk(d1 < 1e-10, sprintf("periodic shift: exported == live (max diff %.2e)", d1),
    sprintf("periodic shift: exported differs from live by %.6f", d1))
d2 <- max(abs(live$registered_curves - expd$registered_curves))
chk(d2 < 1e-10, sprintf("periodic registered curves: exported == live (max diff %.2e)", d2),
    sprintf("periodic registered curves differ by %.6f", d2))

# The regression that was actually shipped: the old export attenuated the shift
# by 0.1 * 0.5 and had no periodic branch. Guard the MAGNITUDE directly, so a
# reintroduced attenuation fails here even if both sides drift together.
got <- abs(live$shifts[3])
chk(abs(got - 0.10) < 0.02,
    sprintf("a true 0.10 displacement is recovered as %.4f", got),
    sprintf("a true 0.10 displacement came back as %.4f", got))

# -- parametric, every family the UI offers -----------------------------------
for (fam in c("power", "exponential", "quadratic", "logistic")) {
  lv <- kern$parametric_alignment(fd, family = fam, param_range = c(0.5, 2), time_points = tp)
  ex <- run_emitted("parametric_alignment", sprintf(
    "warp_result <- parametric_alignment(fd_obj, family = '%s', param_range = c(0.5, 2), time_points = time_points)", fam))
  dd <- max(abs(lv$alpha_values - ex$alpha_values))
  chk(dd < 1e-10, sprintf("parametric/%s: exported alpha == live (max diff %.2e)", fam, dd),
      sprintf("parametric/%s: exported alpha differs by %.6f", fam, dd))
  chk(lv$identity_alpha >= lv$param_range_used[1] && lv$identity_alpha <= lv$param_range_used[2],
      sprintf("parametric/%s: identity alpha = %g is inside the range used", fam, lv$identity_alpha),
      sprintf("parametric/%s: identity alpha %g is OUTSIDE the range used", fam, lv$identity_alpha))
}

# Curves that need no warping must not be deformed — the stale export moved the
# time axis by 0.4999 of the domain here.
Z  <- sapply(1:5, function(i) prof(tp))
fz <- smooth.basis(tp, Z, fdPar(create.fourier.basis(c(0,1), nbasis = 15), 2, 1e-8))$fd
for (fam in c("power", "exponential", "quadratic", "logistic")) {
  r <- kern$parametric_alignment(fz, family = fam, param_range = c(0.5, 2), time_points = tp)
  dev <- max(abs(r$warp_functions - tp))
  chk(dev < 0.02,
      sprintf("parametric/%s leaves identical curves alone (max |h(t)-t| = %.4f)", fam, dev),
      sprintf("parametric/%s deformed curves needing no warp by %.4f of the domain", fam, dev))
}

# -- landmark ------------------------------------------------------------------
lm_live <- kern$landmark_alignment_simple(fd, c(0.25, 0.5, 0.75), tp)
chk(!is.null(lm_live), "landmark registration returns a result",
    "landmark registration returned NULL (the P11.2 scoping defect)")
if (!is.null(lm_live)) {
  lm_exp <- run_emitted("landmark_alignment_simple",
    "warp_result <- landmark_alignment_simple(fd_obj, c(0.25, 0.5, 0.75), time_points, landmark_points = NULL)")
  d3 <- max(abs(lm_live$warp_functions - lm_exp$warp_functions))
  chk(d3 < 1e-10, sprintf("landmark: exported warp == live (max diff %.2e)", d3),
      sprintf("landmark: exported warp differs by %.6f", d3))
  # the helper closure must travel with it, or the script dies on line one
  chk(any(grepl("fck_landmark_warp <- function", emitted)),
      "the exported script carries fck_landmark_warp with its caller",
      "the exported script calls fck_landmark_warp but never defines it")
}

# ---- 5. the generator must not hand-write these algorithms again -------------
gen <- paste(readLines("server/90_export.R", warn = FALSE), collapse = "\n")
# The attenuation constants are the signature of the pre-P0.8 copy. They are
# matched as CODE THE GENERATOR WOULD EMIT, not quoted from the comment that
# explains their removal — a guard that greps for text its own comment contains
# is vacuous, which has happened three times in this audit.
sigs <- c('add\\("    shifts\\[i\\] <- best_lag / n_time \\* 0\\.1"\\)',
          'add\\("    warp_functions\\[,i\\] <- pmin')
for (sig in sigs)
  chk(!grepl(sig, gen), paste("90_export.R does not hand-write:", sig),
      paste("90_export.R still hand-writes the registration:", sig))
for (k in c("linear_shift_alignment", "parametric_alignment", "landmark_alignment_simple"))
  chk(grepl(paste0('emit_kernel\\("', k, '"\\)'), gen),
      paste("90_export.R emits", k, "by deparse()"),
      paste("90_export.R does not emit", k, "from the live definition"))

cat(sprintf("\n%s  (%d passed, %d failed)\n",
            if (bad == 0) "Warp export roundtrip tests PASSED" else "Warp export roundtrip tests FAILED",
            ok, bad))
quit(status = if (bad == 0) 0 else 1)
