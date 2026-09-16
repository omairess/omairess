# ==============================================================================
# tests/traj_tab6_test.R — the Group / Condition Comparison pipeline
# ==============================================================================
# Tab 6 is five panels over ONE fitted model. This drives the same calls the
# reactives make, for each of the three designs, and checks that every panel has
# something to render -- including the panels that must REFUSE.
#
# It also checks the wiring: every output id the UI asks for must exist in the
# server. A missing one renders as a blank space with no error anywhere, which
# is the failure mode a smoke test that only parses files cannot see.
# ==============================================================================

`%||%` <- function(a, b) if (is.null(a)) b else a
app <- if (dir.exists("server")) "." else ".."
e <- new.env()
for (f in c("02b_helpers_palette", "08_helpers_cosinor", "08c_helpers_circstat",
            "08d_helpers_traj", "08e_helpers_trajfit", "08f_helpers_trajinf",
            "08g_helpers_trajband"))
  sys.source(file.path(app, "server", paste0(f, ".R")), envir = e)

pass <- 0L; fail <- 0L
chk <- function(ok, a, b) { if (isTRUE(ok)) { pass <<- pass + 1L; cat("ok   ", a, "\n") }
                            else { fail <<- fail + 1L; cat("FAIL ", b, "\n") } }

# ---------------------------------------------------------------- the wiring
cat("-- every output the UI asks for exists in the server ------------------\n")
ui  <- paste(readLines(file.path(app, "ui/72_harmonic.R"), warn = FALSE), collapse = "\n")
# EVERY server file, not just 72_harmonic.R: the app sources them all into one
# environment, so an output declared in the harmonic UI may legitimately be
# rendered elsewhere -- tab 2b's density outputs live in 74_polar_density.R.
# The first draft of this check searched one file and reported three false
# positives, which is the same class of mistake as the test itself hunts for.
srv <- paste(unlist(lapply(
  list.files(file.path(app, "server"), pattern = "[.]R$", full.names = TRUE),
  readLines, warn = FALSE)), collapse = "\n")
srv_harm <- paste(readLines(file.path(app, "server/72_harmonic.R"), warn = FALSE),
                  collapse = "\n")
ids <- unique(unlist(regmatches(ui, gregexpr(
  '(uiOutput|plotlyOutput|verbatimTextOutput|plotOutput|tableOutput|DTOutput)\\("[A-Za-z0-9_.]+"', ui))))
ids <- sub('^.*\\("', "", sub('"$', "", ids))
missing <- ids[!vapply(ids, function(i)
  grepl(sprintf("output\\$%s\\s*<-", i), srv), logical(1))]
chk(length(missing) == 0,
    sprintf("all %d output ids in the harmonic UI are rendered by the server", length(ids)),
    sprintf("no server output for: %s", paste(missing, collapse = ", ")))

# THE OTHER DIRECTION, and the one that actually bit.
# The check above finds a uiOutput with no renderer. The reverse -- a server
# that reads input$X where X is placed by no UI element -- is worse, because
# Shiny returns NULL for a missing input and any `if (!is.null(...))` guard
# quietly takes the "no" branch. Replacing the Group Variable control with the
# Study Design panel orphaned input$harmonic_group_var, and eleven outputs
# (the fitted-curve plot, the polar plot, four parameter histograms, the
# individual table, the export, the legacy comparison) silently drew ungrouped
# results with no error anywhere.
# EVERY server file, not just the harmonic one. Scoping this to 72_harmonic.R
# let the identical bug survive in 73_cosinor_pairwise.R, which reads the same
# inputs -- the third time in this file that a check was written against too
# narrow a slice of the app. All server files share one environment; the check
# should share its scope.
srv_inputs <- unique(unlist(regmatches(srv, gregexpr('input\\$[A-Za-z0-9_.]+', srv))))
srv_inputs <- sub("^input\\$", "", srv_inputs)
ui_all <- paste(unlist(lapply(
  list.files(file.path(app, "ui"), pattern = "[.]R$", full.names = TRUE),
  readLines, warn = FALSE)), collapse = "\n")
# An input counts as placed only if some file CREATES it -- i.e. its name is the
# first argument of an input constructor. A first draft accepted the name
# appearing in any quoted string, which a comment or a list element satisfies,
# and the check passed against a deliberately re-broken copy. A test that cannot
# fail on the bug it was written for is worse than no test.
# An input counts as placed only if its name is the FIRST ARGUMENT OF A CALL --
# either a shiny constructor directly, or a helper that builds one, since the
# design pickers are rendered by .harmonic_factor_picker(id, label). Two earlier
# drafts got this wrong in opposite directions and were caught by re-breaking a
# scratch copy on purpose: accepting the name in any quoted string passed the
# broken copy (a comment satisfied it), and listing only shiny constructors
# flagged the four dynamically-built pickers as orphans. A check that cannot
# fail on the bug it was written for is worse than no check.
placed <- vapply(srv_inputs, function(i) {
  pat <- sprintf('[A-Za-z0-9_.]+\\(\\s*(inputId\\s*=\\s*)?"%s"', i)
  grepl(pat, ui_all) || grepl(pat, srv)
}, logical(1))
# FOUR PRE-EXISTING ORPHANS, verified benign and recorded rather than hidden.
# These predate P21 and live in other modules. Each read is guarded -- the code
# checks is.null() and takes a defined branch -- so they degrade to a documented
# default rather than silently changing a result, which is what made the
# harmonic_group_var case different. They are listed here so the check can fail
# on a NEW one; removing a name from this list must mean the input was created,
# never that the check became inconvenient.
known_benign <- c(
  min_observed_points = "10_import.R: guarded by is.null/is.finite; no row filter applied when absent",
  selected_subject    = "40_fpca.R: guarded; falls back to the mean landmark target",
  landmark_target     = "40_fpca.R: guarded; same fallback",
  density_avg_days    = "74_polar_density.R: legacy fallback, only read when density_profile_mode is absent")
orphans <- setdiff(srv_inputs[!placed], names(known_benign))
chk(length(orphans) == 0,
    sprintf("all %d inputs the server reads are placed by some UI", length(srv_inputs)),
    sprintf("the server reads inputs nothing creates: %s", paste(orphans, collapse = ", ")))

traj_ids <- grep("^harmonic_traj", ids, value = TRUE)
chk(length(traj_ids) >= 8,
    sprintf("tab 6 contributes %d of them (%s...)", length(traj_ids),
            paste(utils::head(traj_ids, 3), collapse = ", ")),
    "tab 6's outputs are missing from the UI")

# ---------------------------------------------------------------- fixtures
gen <- function(seed, groups = c("a","b"), conds = NULL, n_per = 6, nt = 12,
                slope = 0, amp = 4) {
  set.seed(seed); tp <- seq(0, 22, length.out = nt)
  cs <- conds %||% NA_character_
  subj <- sprintf("S%03d", seq_len(n_per * length(groups)))
  sg <- rep(groups, each = n_per)
  rows <- list(); k <- 1L
  meta <- list(subject = character(0), Group = character(0), Condition = character(0))
  for (i in seq_along(subj)) {
    sb <- rnorm(1, 0, 2); sa <- rnorm(1, 0, .4)
    for (cc in cs) {
      rows[[k]] <- 10 + sb + rnorm(1,0,1.5) + (amp + sa) * cos(2*pi*tp/24 - 2) +
                   slope * tp + rnorm(nt, 0, 1)
      meta$subject <- c(meta$subject, subj[i]); meta$Group <- c(meta$Group, sg[i])
      meta$Condition <- c(meta$Condition, cc); k <- k + 1L
    }
  }
  list(Y = do.call(rbind, rows), t = tp, meta = meta)
}
build <- function(g, within) {
  fl <- list(Group = g$meta$Group)
  if (within) fl$Condition <- g$meta$Condition
  e$dance_traj_long(g$Y, g$t, g$meta$subject, fl)
}

# Everything tab 6 renders, for one fit. Mirrors the reactives one to one.
drive <- function(ff, label) {
  cat(sprintf("\n-- %s\n", label))
  chk(isTRUE(ff$ok), sprintf("the model fits: %s", ff$re_formula),
      sprintf("the fit failed: %s", ff$message %||% "?"))
  if (!isTRUE(ff$ok)) return(invisible())

  # 1. primary -- one full-trajectory test per factorial effect
  dts <- ff$spec$design_terms
  tl <- attr(stats::terms(stats::as.formula(ff$spec$fixed_formula)), "term.labels")
  effects <- unlist(lapply(seq_along(dts), function(k)
    utils::combn(dts, k, paste, collapse = ":")), use.names = FALSE)
  got <- vapply(effects, function(ef) {
    parts <- strsplit(ef, ":", fixed = TRUE)[[1]]
    keep <- tl[vapply(tl, function(tm) setequal(
      intersect(strsplit(tm, ":", fixed = TRUE)[[1]], dts), parts), logical(1))]
    r <- e$dance_traj_block_test(ff, keep)
    isTRUE(r$ok) && is.finite(r$p)
  }, logical(1))
  chk(all(got), sprintf("a full-trajectory test for each of %d effect(s): %s",
                        length(effects), paste(effects, collapse = ", ")),
      sprintf("no testable block for: %s", paste(effects[!got], collapse = ", ")))

  # 2. secondary -- the component decomposition
  set <- e$dance_traj_omnibus_set(ff, which = c("full","shape","circadian","trend","level"))
  okb <- vapply(set$results, function(r) isTRUE(r$ok), logical(1))
  expect_trend <- !identical(ff$spec$trend, "none")
  chk(all(okb[c("full","shape","circadian","level")]) &&
        identical(unname(okb["trend"]), expect_trend),
      sprintf("the decomposition returns %d blocks; the trend row appears only with a trend (%s)",
              sum(okb), ff$spec$trend),
      sprintf("block availability is wrong: %s", paste(names(okb)[!okb], collapse = ", ")))

  # WITHOUT A TREND, two of the five blocks are the same test and two of the
  # four component views are degenerate. Found by running on real data, where
  # shape and circadian both came back F(12, 1823.5) = 2.028 -- one test printed
  # twice under two names reads as corroboration. Tab 6 now shows only the
  # blocks and views that differ; this pins the fact that justifies it.
  if (identical(ff$spec$trend, "none")) {
    chk(identical(set$results$shape$terms, set$results$circadian$terms) &&
          abs(set$results$shape$statistic - set$results$circadian$statistic) < 1e-12,
        "with no trend, shape and circadian ARE the same test (so only one is shown)",
        "shape and circadian differ with no trend in the model")
    # It used to return a curve of exact zeros. It now REFUSES: a view that can
    # only draw a flat line is not a component of this model, and offering it
    # and then drawing nothing is how two of the four selector entries came to
    # be a zero line and a duplicate.
    tv0 <- e$dance_traj_predict(ff, component = "trend", n_time = 10)
    chk(!isTRUE(tv0$ok) && grepl("not a component of this fit", tv0$message %||% ""),
        "and the nonperiodic view is refused outright, not drawn as zeros",
        "the nonperiodic view was accepted on a model with no trend")
    chk(!any(c("trend", "baseline_harm") %in% e$dance_traj_components(ff)),
        "neither trend view is offered for this fit",
        "a trend view is still offered with no trend fitted")
  }

  # 3. the component views the FIT offers, each with a band
  comps <- e$dance_traj_components(ff)
  pv <- lapply(comps, function(cc) e$dance_traj_predict(ff, component = cc, n_time = 40))
  names(pv) <- comps
  chk(all(vapply(pv, function(p) isTRUE(p$ok) && all(is.finite(p$table$fit)), logical(1))),
      sprintf("all %d component views render, per design cell", length(comps)),
      "a component view failed")
  # the identity that makes them views of ONE model rather than several models
  if (all(c("baseline_harm", "trend") %in% comps)) {
    fl <- pv$full$table; bh <- pv$baseline_harm$table; tr <- pv$trend$table
    chk(max(abs(fl$fit - (bh$fit + tr$fit))) < 1e-8,
        sprintf("full == baseline+harmonics plus trend (max diff %.1e)",
                max(abs(fl$fit - (bh$fit + tr$fit)))),
        "the component views do not decompose the full trajectory")
    t0 <- tr$t == min(tr$t)
    chk(max(abs(tr$fit[t0])) < 1e-9 && max(abs(tr$hi[t0] - tr$lo[t0])) < 1e-9,
        "the nonperiodic view is zero at the reference time, band included",
        "the trend view is not pinned to zero at the origin")
  }
  # and each harmonic on its own sums back to the harmonics view
  hk <- grep("^harmonic[0-9]+$", comps, value = TRUE)
  if (length(hk)) {
    tot <- Reduce(`+`, lapply(hk, function(k) pv[[k]]$table$fit))
    chk(max(abs(tot - pv$harmonics$table$fit)) < 1e-9,
        sprintf("the %d per-harmonic views sum to 'harmonics only'", length(hk)),
        "the per-harmonic views do not sum to the harmonics view")
  }

  # 4. difference curve
  cells <- e$dance_traj_cell_grid(ff$spec)$.cell
  if (length(cells) > 1) {
    dc <- e$dance_traj_diff_curve(ff, cells[1], cells[2], n_time = 40)
    chk(isTRUE(dc$ok) && all(is.finite(dc$table$diff)),
        sprintf("the difference curve %s - %s renders", cells[1], cells[2]),
        "the difference curve failed")
  }

  # 5. pairwise, and the simple effects the selector offers
  for (w in c("level", "amplitude", "phase")) {
    r <- e$dance_traj_contrasts(ff, w, adjust = "holm")
    chk(isTRUE(r$ok) && nrow(r$table) == choose(length(cells), 2),
        sprintf("pairwise %s: %d contrasts over %d cells", w, nrow(r$table %||% data.frame()), length(cells)),
        sprintf("pairwise %s failed", w))
  }
  if (length(dts) > 1) {
    others <- setdiff(dts, dts[1])
    at <- as.list(stats::setNames(vapply(others, function(f)
      levels(ff$spec$data[[f]])[1], character(1)), others))
    se <- e$dance_traj_simple_effects(ff, dts[1], at = at, what = "level")
    chk(isTRUE(se$ok) && nrow(se$table) >= 1,
        sprintf("simple effect of %s at %s: %d contrast(s)", dts[1],
                paste(sprintf("%s=%s", names(at), unlist(at)), collapse=","), nrow(se$table)),
        "the simple-effect slice failed")
  }
}

suppressWarnings(suppressMessages({
drive(e$dance_traj_fit(e$dance_traj_spec(build(gen(1), FALSE), 24, 1, "none")),
      "BETWEEN: 2 groups, 12 participants")
drive(e$dance_traj_fit(e$dance_traj_spec(
        build(gen(2, groups = "a", conds = c("p","q","r"), n_per = 12), TRUE), 24, 1, "none")),
      "WITHIN: 3 conditions, 12 participants")
drive(e$dance_traj_fit(e$dance_traj_spec(
        build(gen(3, groups = c("a","b"), conds = c("p","q"), slope = .1), TRUE), 24, 1, "linear")),
      "MIXED: 2 x 2 with a linear trend")
drive(e$dance_traj_fit(e$dance_traj_spec(
        build(gen(4, groups = c("a","b","c"), conds = c("p","q"), n_per = 5), TRUE), 24, 2, "none")),
      "MIXED: 3 x 2, two harmonics")
}))

cat(sprintf("\n%s  (%d passed, %d failed)\n",
            if (fail == 0) "Tab 6 pipeline tests PASSED" else "FAILURES", pass, fail))
if (fail > 0) quit(status = 1)
