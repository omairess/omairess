# WINSTEPPER — project memory

A Shiny app for WINSTEPS-style Rasch measurement. **This is a fresh, redesigned
front end built on the same estimation engine as R-Winsteps.** The estimation
engine is written from scratch in base R so every number traces to a published
formula. It is self-contained in this folder and does not touch the sibling
BootSON / Dagger / PsychoNetrix work in the repo.

## Layout

| File | Role |
|---|---|
| `rasch_engine.R` | JMLE estimation, fit statistics, category structure, DIF, PCA of residuals. **No Shiny dependency** — must stay usable from the command line. **Reused unchanged from R-Winsteps.** |
| `winsteps_plots.R` | All figures, base graphics, every one driven by a `style` list. **Reused unchanged.** |
| `house_modules.R` | Shared house-style modules (data load, reshape, varselect, step recorder, exports). **Reused unchanged.** |
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

Everything below the UI — the engine API, the recorder discipline, the input IDs
the server binds to — is deliberately identical to R-Winsteps so the numbers and
the reproducible export stay the audited ones.

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

**Model-expected correlations are not `cor(E, θ)`.** Under the model
`Var(X) = Var(E) + mean(W)`, so the expected point-measure correlation is
`cov(E, θ) / sqrt((Var(E) + mean(W)) * Var(θ))`.

## Statistical scope — be honest in any output

- Not byte-identical to WINSTEPS 5.11. Agreement ≈ 2 decimals on well-behaved data.
- No `STBIAS=` JMLE bias correction, so the usual spread inflation (~L/(L−1)) is
  present. WINSTEPS shows it too by default.
- Per-category INFIT MNSQ is not centred on 1.0 in extreme categories. Linacre's
  "< 2.0" guideline applies to the OUTFIT column.
- PCA clusters are tertiles of first-contrast loadings — an approximation to the
  WINSTEPS three-cluster split.
- Not implemented: anchoring (`IAFILE=`/`PAFILE=`/`SAFILE=`), `CUTLO=`/`CUTHI=`,
  subset/connectivity detection, keyform and scalogram tables, DPF (Table 31),
  non-uniform DIF, multi-facet models.
- Independent software; not affiliated with or derived from WINSTEPS®.

Never loosen a statistical claim to make a test pass. If a result cannot be
reproduced by the exported script, that is a blocking bug in the recorder wiring.

## Likely next work

1. Anchoring (`IAFILE=`/`PAFILE=`/`SAFILE=`) — the biggest gap; needed for equating
   across forms or waves.
2. `STBIAS=`-style bias correction as an opt-in.
3. DPF (Table 31) and non-uniform DIF.
4. Subset/connectivity detection — currently a silent failure mode on sparse data.
5. Port the R-Winsteps test suites into this folder and wire them to `app.R`.
