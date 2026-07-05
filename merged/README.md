# Merged network-psychometrics workbench

One Shiny app combining three previously-standalone tools:

- **GGM (bootnet)** — regularized/non-regularized Gaussian graphical models, from `BootSON.R`
- **DAG (bnlearn)** — Bayesian network structure learning, from `Dagger_zero.R`
- **Latent (psychonetrics) & NCT** — GGM/LNM/RNM/LRNM models and network comparison, from `PsychoNetrix.R`

Built to `shiny-house-style` (ten house rules, step-recorder pattern) and
`network-psychometrics-correctness` (mandatory bootstraps, CS-coefficient
gate, CPDAG caveat, paired NCT + multiple-comparison correction). See
`stage1_inventory_report.md` (repo root) for the source-app audit and
`stage2_decisions.md` for the integration decisions.

## Layout

```
merged/
  app.R                        — assembly: navbar UI + server wiring
  R/
    00_packages.R               shared  — package guard + version ledger (rule 0)
    01_recorder.R               shared  — THE step recorder (rules 6 + 8)
    02_palette.R                shared  — pastel default colours (rule 9)
    03_data_module.R            shared  — file loader + wide/long reshape (rules 1, 2)
    04_varselect_module.R       shared  — searchable multi-select variable picker (rule 3)
    05_appearance_module.R      shared  — colour pickers, size sliders, PNG/PDF/SVG export (rules 4, 5)
    10_tab_bootnet.R            per-tab — bootnet GGM estimation + mandatory bootstraps + CS gate
    11_tab_dagger.R             per-tab — bnlearn structure learning + boot.strength + CPDAG caveat
    12_tab_psynet.R             per-tab — psychonetrics models + NCT (paired/independent, corrected)
```

### Shared modules (used by all three tabs)

Every tab consumes the **same** data contract and the **same** recorder
instance — there is one canonical dataset and one exported script/pipeline
log for the whole app, not one per tab.

- **Data contract** (`03_data_module.R`): `data_bus$raw()` is the file exactly
  as uploaded; `data_bus$wide()` is the canonical analysis frame every tab
  reads from (long data gets pivoted to wide with `var_<time>` column
  naming — chosen because the DAG tab's temporal-pattern detection greps for
  that suffix). Both raw and pivoted frames get a head+dims preview.
- **Step recorder** (`01_recorder.R`): every analysis stage in every tab
  registers `{id, phase, description, code}`. Re-running a stage upserts by
  id. The "Pipeline & Export" tab's plain-language log and downloadable R
  script are both rendered from this one list, so they can never drift out
  of sync with each other or with what the app actually did.
- **Variable picker** (`04_varselect_module.R`), **appearance module**
  (`05_appearance_module.R`: colour pickers, size sliders, PNG/PDF/SVG
  export), **palette** (`02_palette.R`: RColorBrewer Pastel1-derived
  defaults), **package ledger** (`00_packages.R`: install-missing-only,
  version-recording) are each a single implementation all three tabs share,
  replacing three divergent ad hoc versions in the source apps.

### Per-tab modules (analysis-family-specific, not shared)

Each tab owns its own model-fitting logic, its own stability/validation
method (these are genuinely different statistical procedures, not
interchangeable), and its own results display:

- **GGM tab**: estimator dispatch (EBICglasso primary), nonparametric +
  case-dropping bootstraps, CS-coefficient gate, expected-influence-vs-strength
  logic.
- **DAG tab**: bnlearn algorithm/score dispatch, `boot.strength` +
  `averaged.network`, CPDAG/Markov-equivalence audit and caveat.
- **Psychonetrics/NCT tab**: model-family dispatch, estimator resolution
  (ordinal/ML mismatch guard), NCT with explicit paired handling and
  corrected edge-level tests (pinned to the GGM tab's own estimator settings,
  so a comparison can never silently use different settings than the
  networks on screen).

## What's implemented vs. still a stub

Stage 2 wrote the correctness-critical logic above in full (estimator
pinning, mandatory bootstraps, the CS gate, the CPDAG caveat, paired NCT +
correction). Stage 3 filled in the generic mechanical scaffolding: the
extension-dispatch file loader (`.csv/.txt` with delimiter + decimal-mark
control, `.tsv`, `.xls/.xlsx` with sheet selection, `.sav`, `.rds`,
`.RData` with an object picker), the searchable variable picker, the
colour-picker + slider appearance controls wired into the GGM and DAG tabs'
plots (with PNG/PDF/SVG export re-plotting into a fresh device), and the
pastel default palette.

Added in the 2026-07-05 change requests:

- GGM tab: on-demand bootstrap validation (instant estimation; a separate
  button runs the bootstraps — the CS-coefficient gate still locks
  centrality until it has run for the current network); all eight
  estimators (EBICglasso / ggmModSelect / pcor / cor / huge / IsingFit /
  mgm / relimp); layout options (spring / circle / EGA communities / PCA);
  node-size-by-column-mean.
- DAG tab: blacklist/whitelist constraint builder (threaded into every
  bootstrap replicate); advanced arc metrics (strength / direction /
  combined) driving width and labels; layout options (spring / circle /
  hierarchical Sugiyama); node-size-by-column-mean.

Still open — **larger feature ports from the source apps**, each marked
`TODO(port)` at its point of use:

- GGM tab: Bridge-Symptoms centrality sub-tab, the NSON sub-app,
  split-group estimation, bootstrap edge-significance thresholding.
- DAG tab: split analysis, DAG-diff plot, the folded-temporal-graph
  feature, the full 26-entry bnlearn score list.
- Psychonetrics tab: non-GGM model families (LNM / RNM / LRNM / panel
  models), the lambda-matrix editor, prune/step-up/model-search controls,
  and the **latent** Relative-Importance network (Johnson/LMG) — distinct
  from the GGM tab's observed-variable `relimp` estimator.
- NCT: subject-ID matching UI for paired designs (currently requires
  equal-sized, row-aligned groups as an interim guard).

None of these are wired into the shared recorder yet, so none of them can
silently violate either governing skill in the meantime — they simply
aren't present.
