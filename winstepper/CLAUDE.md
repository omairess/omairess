# WINSTEPPER — project memory

A Shiny app for WINSTEPS-style Rasch measurement. **This is a fresh, redesigned
front end built on the same estimation engine as R-Winsteps.** The estimation
engine is written from scratch in base R so every number traces to a published
formula. It is self-contained in this folder and does not touch the sibling
BootSON / Dagger / PsychoNetrix work in the repo.

## Layout

| File | Role |
|---|---|
| `rasch_engine.R` | JMLE estimation, fit statistics, category structure, DIF, PCA of residuals. **No Shiny dependency** — must stay usable from the command line. **Byte-identical to R-Winsteps; do not edit it.** |
| `winstepper_cmle.R` | Conditional maximum likelihood (CMLE): elementary symmetric functions, conditional likelihood + analytic gradient, `rasch_cmle()`. **No Shiny dependency.** WINSTEPPER-only; source it **after** the engine. |
| `winsteps_plots.R` | All figures, base graphics, every one driven by a `style` list. Originally from R-Winsteps; now also holds the label-fitting helpers. |
| `house_modules.R` | Shared house-style modules (data load, reshape, varselect, step recorder, exports). **Reused unchanged.** |
| `winstepper_extras.R` | WINSTEPPER-only additions kept out of the audited files: general keyform (Table 2.2) and DGF (Table 33). Depends on engine + plot internals; source it **after** them. |
| `app.R` | The WINSTEPPER Shiny app — a redesigned `bslib` UI (grouped nav menus, value-box dashboard, cards) over the engine above. |

## What is "fresh" about WINSTEPPER vs R-Winsteps

- UI is `bslib::page_navbar` with a Bootstrap 5 theme, not the old flat
  `navbarPage`. The 13 tabs are grouped under **Data**, **Estimate**,
  **Results** (Summary / Wright map / Items / Persons / Rating scale / Score
  table), **Advanced** (Graphs / Dimensionality / DIF) and a right-aligned
  **Settings** menu (Figure style / Exports / About).
- A value-box dashboard surfaces persons/items/iterations/convergence on the
  Estimate tab and person/item reliability + separation on the Summary tab.
- Panels use `bslib::card` / `layout_columns` / `layout_sidebar`.

## WINSTEPPER-only features (not in R-Winsteps)

- **CODES= / NEWSCORE=** on the Estimate tab. `CODES=` declares which response
  codes are valid categories (others → missing); `NEWSCORE=` positionally
  recodes them (e.g. `0,1,2,3` → `0,1,1,2` collapses categories 1 and 2).
  Applied before recoding; both are written into the reproducible script
  (`apply_codes` / `parse_newscore` / `apply_codes_recode` in app.R).
- **Grouped / mixed model scale editor.** Choosing the grouped model shows a
  per-scale editor: assign items to each scale via a `pickerInput` and give each
  scale its own `CODES` / `NEWSCORE` (blank = use the global defaults). This
  builds the `groups` vector passed to `rasch_prep()`.
- **Category-count-aware threshold advance.** `threshold_advance_min(m)` in
  `winstepper_extras.R` returns `ln((k+1)(m+1-k) / (k(m-k)))` for k = 1..m-1.
  The familiar "advance by 1.4 logits" is only the **three-category** case
  (2 ln 2 = 1.386); with more categories the requirement is smaller (4 cats:
  1.10, 1.10; 5 cats: 0.98, 0.81, 0.98; …). Verified against Linacre's table
  (RMT 2006, 20:1, p. 1052) for 3–10 categories.
  `winstepper_extras.R` **deliberately overrides** the engine's
  `category_diagnostics()` to use this (the engine tested a flat 1.4). Since
  extras is sourced last, the override wins — do not "fix" this by editing
  `rasch_engine.R`, which is kept byte-identical to R-Winsteps.
- **General keyform (Table 2.2)** — `keyform_data()` / `plot_keyform()` in
  `winstepper_extras.R`, with expected-score (2.2), Rasch-Thurstone (2.3) and
  modal (2.1) variants, and an optional **person-measure histogram** drawn
  underneath on the *same* logit axis (`show_persons`). Records a `keyform` step.
- **Person include/exclude + re-estimation.** `input$excl_persons` holds the
  exclusion list; `prep()` drops those rows and stores `keep_rows` / `excluded`
  on the prep object. **Anything that pairs a column of `wide_data()` with the
  fit must subset by `prep()$keep_rows`** — DIF and DGF do; add the same to any
  new person-level covariate analysis or the lengths will silently mismatch.
- **Interface font size.** `input$ui_font` / `input$table_font` render into
  `output$font_css`, which sets the root `html{font-size}` (Bootstrap is
  rem-based, so the whole UI scales together).
- **Person-item barchart.** `pi_threshold_data()` / `plot_pi_map()` in
  `winstepper_extras.R` — the WINSTEPS Plots-menu chart: person-measure
  histograms (whole sample plus one panel per level of a classification) beside
  every item drawn as a vertical string of numbered thresholds, all on one
  shared logit axis. Exposes the WINSTEPS options: measure-relative values
  (Andrich / Thurstonian / half-point / maximum-probability full-point), sort
  order (entry / measure / alpha, ascending or descending), right-side size and
  title, and "add item measures". Lives on the Wright map tab; records a `pimap`
  step. **Its `person_class` must be subset by `prep()$keep_rows`**, like DIF and
  DGF, or the lengths silently disagree. All panels must keep identical top and
  bottom margins or the boxes stop lining up under `layout()`.
- **Suggested rescore.** `suggest_collapse()` in `winstepper_extras.R` proposes a
  `NEWSCORE=` collapse for one item group. **It is a suggestion only**: the button
  on the Rating scale tab fills the boxes and navigates to the Estimate tab, but
  never re-estimates — the user confirms or edits and presses Estimate.
  Deliberately no LLM/API: it must be free, offline, deterministic and
  reproducible. Covered by `test_collapse.R`.
  **When filling the boxes, always write CODES too** — `apply_codes_recode()`
  ignores NEWSCORE when CODES is blank, so writing NEWSCORE alone silently does
  nothing.

  **Do not re-introduce a "merge until the advance guideline is satisfied" rule.**
  The first version did that and it failed badly on real data
  (`NSGGMgameCONT.xlsx`, 17 items × 1025 persons, 10 categories, none sparse):
  when the thresholds are genuinely compressed no collapse can ever satisfy the
  guideline, so the rule ran to the `min_cat` floor and proposed `0 | 1-3 | 4-9`
  — two poles and one enormous middle, which guts the information function and
  person separation. The objective is **inverted**: keep as MANY categories as
  possible subject to their thresholds staying `min_sep` logits apart, tie-broken
  by even spacing. On that dataset the corrected version returns
  `0 | 1-3 | 4-6 | 7-8 | 9`, matching independent expert judgement.

  The basis is `cum_threshold_locations()`: `-logit(P(X >= k))` from the observed
  counts. These have an exact subset property under adjacent merging — the
  thresholds of any candidate collapse are literally a subset of these — so every
  candidate can be scored with **no re-estimation and no approximation**. The
  search is exhaustive (2^m partitions, 1024 for 11 categories), with a greedy
  fallback above `max_enum_cat`. `min_sep` is exposed in the UI: raise it for a
  coarser scale, lower it to retain more categories.
- **CMLE** — `rasch_cmle()` in `winstepper_cmle.R`, selected by `input$method`
  on the Estimate tab (JMLE stays the default). Conditions each response pattern
  on its person raw score, so the person parameters leave the likelihood and the
  item estimates are consistent — no `STBIAS=`-style correction needed. Returns
  the **same `raschfit` structure** `rasch_jmle()` does, plus a `cml` sub-list
  and `settings$method`, so every downstream table and figure works untouched.
  Points that are easy to get wrong, and are the reason `test_cmle.R` exists:
  - The elementary symmetric functions **must** be computed in the log domain.
    On the NSGGM data (17 items × 10 categories, raw scores to 153) the raw
    gammas span hundreds of orders of magnitude and overflow instantly.
  - `G^(-i)` is obtained by **re-running the convolution without item i**, never
    by the difference/deconvolution algorithm — that is faster but numerically
    unstable, and at these test lengths the safe route is free.
  - Missing data conditions **within each missingness pattern**: a person's
    conditioning set is the items they actually answered, so each pattern needs
    its own gamma. This is the main cost driver on sparse data.
  - Parameterised as the engine is (`b_ik = k·delta_i + cumsum(tau_g)[k]`) with
    sum-to-zero imposed *by construction*, not by projection, so BFGS has no
    flat direction to wander along.
  - BFGS stops on the function value and leaves the score residual near 1e-4.
    Convergence is judged on the CML equations (observed = conditionally
    expected sufficient statistic) and finished with Newton steps using the
    Hessian that the standard errors need anyway.
  - **Read `fit$settings$method`, never `input$method`**, anywhere that describes
    a result — the radio can be changed without pressing Estimate.
- **DGF (Table 33)** — `dgf_analysis()` / `plot_dgf()`. Item classes come from
  the model scales (`fit$groups`); one uniform difficulty shift is estimated per
  item-class × person-class cell and contrasted across person classes with the
  ETS A/B/C rule. Records a `dgf` step. Sits alongside DIF under **Advanced**.

New step ids `keyform` and `dgf` are already in `STEP_ORDER`.

**Estimation is driven by a `refit` counter**, not by `input$run` directly: both
"Estimate" and "Re-estimate without excluded persons" bump it, and `prep()` /
`fit()` are `eventReactive(refit(), ignoreInit = TRUE, ...)`. `ignoreInit` is
required — a plain `reactiveVal(0)` is *not* a "null event" (only `NULL` and an
actionButton at 0 are), so without it both would evaluate at startup.

Everything below the UI — the engine API, the recorder discipline, the input IDs
the server binds to — is deliberately identical to R-Winsteps so the numbers and
the reproducible export stay the audited ones.

## The "reused unchanged" rule, narrowed

`rasch_engine.R` stays **byte-identical** to R-Winsteps. That is what protects
the audited numbers, and it is why `winstepper_extras.R` *overrides* the
engine's `category_diagnostics()` / `category_table()` / `.summary_block()`
rather than fixing them in place.

`winsteps_plots.R` is no longer under that rule. It contains no estimation
logic, and the label-fitting work had to touch four functions in it; copying
them into extras as overrides would have duplicated ~160 lines that then have to
be kept in sync. Changes there must still be justified and covered by
`test_collapse.R`.

## Running

```r
shiny::runApp("winstepper")   # from the repo root
# or, from inside winstepper/:
shiny::runApp("app.R")
```

Packages: shiny, **bslib**, shinyWidgets, colourpicker, DT, readr, readxl,
haven, tidyr, tidyselect, RColorBrewer. `bslib` is the only addition over
R-Winsteps.

> The engine/plots/modules are copied verbatim from the audited R-Winsteps
> project. If you port the R-Winsteps test suites (`tests_engine.R`,
> `test_app.R`, `test_plots.R`) into this folder, run all three before declaring
> any change done — every engine bug found so far was found by them, not by
> reading. `test_app.R` will need its input IDs checked against this app.R (the
> IDs are the same, but the tab structure differs).

## Non-negotiable conventions (house style)

0. Ensure packages, record versions in the export. Never auto-update.
1. Accept `.txt/.csv/.xls/.xlsx/.sav/.rds/.RData`; always preview head + dims first.
2. Wide↔long conversion is an explicit, inspectable step.
3. Variable selection is a searchable multi-select (`pickerInput`).
4. Every figure exposes colour pickers.
5. Every figure exposes size sliders and PNG/PDF/SVG export.
6. **The app always emits a standalone R script that reproduces the analysis.**
7. Always allow export of data and result objects.
8. Always show an ordered, plain-language pipeline log.
9. Pastel palettes by default (`RColorBrewer` Pastel1/Pastel2).

Rules 6 and 8 share one mechanism — the step recorder. A stage may compute a
result **only if** it also records the code that reproduces it. No orphan
computations.

## Traps that have already bitten — do not reintroduce

**`NA * FALSE` is `NA` in R — masking by multiplication does not zero out NAs.**
`rasch_engine.R:683` does `w <- P[[k+1]] * selg; if (sum(w) > 0)`. `P` contains
NA for persons with no estimable measure, and multiplying by the logical mask
leaves those NAs in place, so `sum(w)` is NA and the `if` aborts. Mask by
**assignment** (`w[!selg] <- 0`), never by multiplication, whenever the operand
can contain NA.

**Zone bounds can be NA.** `rasch_engine.R:705` does `if (sum(inzone) > 0)` where
`inzone` depends on `.score_to_measure()`, which returns NA when the root is not
bracketed in ±25 logits. `TRUE & NA` is NA, so `inzone` is poisoned. Long rating
scales (0–10) hit this because the thresholds spread far enough that the
half-score points fall outside the bracket. `winstepper_extras.R` overrides
`category_table()` (and `threshold_data()`) with the mask-by-assignment fix, an
NA-zone guard, and a widened ±60 bracket via `.sm_wide()` / `.th_wide()` —
widening is safe because the expected score is strictly monotone, so the root is
unique and only NAs can change.

**`if (<NA>)` in the engine's summary block.** `rasch_engine.R:579` does
`sep_m <- if (rmse_m > 0) ...`, which throws *"missing value where TRUE/FALSE
needed"* the moment any measure or S.E. in the block is `NA` — killing every
output on the Summary tab at once. This is not rare: `rasch_jmle()` runs its
extreme-**person** loop before the extreme-**item** loop, so a person whose
observed responses all fall on extreme items is skipped (`if (!any(obs)) next`)
and keeps `theta = NA`; the "(all)" rows of Table 3.1 include those cases. Long
rating scales (0–10) trigger it often because floor/ceiling items are common.
`winstepper_extras.R` overrides `.summary_block()` with an NA-tolerant version
that summarises the estimable cases and reports `N_Not_Estimable`;
`measure_health()` + the Summary tab's "Estimability check" card explain who was
dropped. Prefer `isTRUE(x > 0)` over `x > 0` in any new scalar `if`.

**The recorder is not reentrant.** `new_recorder()$record()` reads the `steps`
reactiveVal before writing it. Calling it from a plain `observe()`, `reactive()`
or any `render*()` body creates a reactive dependency on the value it writes →
infinite self-invalidating loop → memory exhaustion. Record **only** from
isolated contexts: `observeEvent()` handlers or `eventReactive()` bodies.
`record()` in `house_modules.R` wraps the read in `shiny::isolate()`; keep it
that way.

**Observer firing order is not guaranteed.** Recorded steps must be re-sorted
into `STEP_ORDER` (defined at the top of the server function) or the exported
script can plot before it estimates. Add any new step id to `STEP_ORDER`.

**The house `mod_varselect_server` emits `data <- data[, vars]`**, which drops
label and DIF columns. It is deliberately passed `null_rec` here; the `prep`
step records the item selection non-destructively instead.

**Wright map panels must share one logit axis.** Compute the histogram breaks
over the *combined* person+item range, not from `hist()` defaults.

**Never hard-code a margin width.** `plot_keyform()` used `par(mar = c(., 8, ., .))`
and clipped `GamingObsession_prob` down to `gObsession_prob`; raising the Figure
style label size made every such figure worse. Size margins with `.fit_labels()`
and pad in-region labels with `.pad_range()` (both in `winsteps_plots.R`).
Two traps in doing so:
- `strwidth()` needs an established plot region, but `mar` must be set *before*
  `plot()`. Hence the `par("cin")`-based estimate, which is deliberately biased
  high — over-estimating wastes margin, under-estimating clips. `par("cin")[1]`
  is **not** the widest glyph (devices set it to 0.9 × pointsize), so weights
  above 1.0 are correct for `W`/`@`/`m`, not a bug.
- Under `layout()`, `par("fin")` and `par("pin")` still describe the *previous*
  panel until the next frame starts. `plot_wright()` therefore passes `fig=` and
  `pin=` computed from the layout widths it set itself. Anything new under
  `layout()` must do the same.

**Model-expected correlations are not `cor(E, θ)`.** Under the model
`Var(X) = Var(E) + mean(W)`, so the expected point-measure correlation is
`cov(E, θ) / sqrt((Var(E) + mean(W)) * Var(θ))`.

## Statistical scope — be honest in any output

- Not byte-identical to WINSTEPS 5.11. Agreement ≈ 2 decimals on well-behaved data.
- No `STBIAS=` JMLE bias correction, so the usual spread inflation (~L/(L−1)) is
  present in a JMLE fit. WINSTEPS shows it too by default. **CMLE** does not have
  it — but our CML is an independent implementation and will not match WINSTEPS'
  CMLE to the last decimal either. Its item S.E.s come from the observed
  information of the *conditional* likelihood, so they do not carry
  person-estimation uncertainty.
- Per-category INFIT MNSQ is not centred on 1.0 in extreme categories. Linacre's
  "< 2.0" guideline applies to the OUTFIT column.
- PCA clusters are tertiles of first-contrast loadings — an approximation to the
  WINSTEPS three-cluster split.
- DGF (Table 33) is a Rasch-Welch-style summary (one uniform shift per cell),
  not WINSTEPS' exact log-linear DGF; treat it as an effect-size indicator.
- Not implemented: anchoring (`IAFILE=`/`PAFILE=`/`SAFILE=`), `CUTLO=`/`CUTHI=`,
  subset/connectivity detection, scalogram table, DPF (Table 31),
  non-uniform DIF, multi-facet models.
- Independent software; not affiliated with or derived from WINSTEPS®.

Never loosen a statistical claim to make a test pass. If a result cannot be
reproduced by the exported script, that is a blocking bug in the recorder wiring.

## Tests

```r
Rscript test_collapse.R   # rescore suggestions, barchart, label fitting
Rscript test_cmle.R       # conditional maximum likelihood
```
Run **both** before declaring any change done. `test_cmle.R` sections 1–3 are
exact checks against brute-force enumeration, a closed form, and a numerically
differentiated likelihood; if CML is edited, those are the ones that catch it.

## Likely next work

1. Anchoring (`IAFILE=`/`PAFILE=`/`SAFILE=`) — the biggest gap; needed for equating
   across forms or waves.
2. `STBIAS=`-style bias correction as an opt-in for JMLE (CMLE now covers the
   case where the bias is what matters).
3. DPF (Table 31) and non-uniform DIF.
4. Subset/connectivity detection — currently a silent failure mode on sparse data.
5. Port the R-Winsteps test suites into this folder and wire them to `app.R`.
