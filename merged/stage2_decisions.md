# Stage 2 — integration decisions (one-line audit notes)

Each non-obvious choice, for review. File references point at the `DECISION:` comments in code.

## Architecture

1. **Three Shiny modules, not renamed IDs** — `NS()` namespacing structurally kills the `run`/`nodeBorderColor`/`layout` collisions and lets each tab keep its own layout vocabulary (semantic conflict dissolved, not papered over). (`app.R`)
2. **One app-wide recorder instance** — rules 6+8 demand a single exported script/log describing everything the app did; per-tab recorders would fragment it. (`01_recorder.R`)
3. **Script fragments ordered phase-first, insertion-second** — guarantees the export runs top-to-bottom in a clean session even when the user hops between tabs out of order. (`01_recorder.R`, unit-tested)
4. **`set.seed(<seed>)` re-emitted before every stochastic fragment**, not once at the top — reproduces the app's numbers regardless of which stages a colleague re-runs. (`01_recorder.R`)
5. **Tab fragments bake their variable subset as an explicit `c("V1",...)` literal** — editing the shared selection after running tab A can't silently redefine tab A's exported code. (`01_recorder.R`, `04_varselect_module.R`)
6. **New data upload drops every downstream recorder fragment (all tabs)** — an exported script must never mix fragments computed from different datasets. (`03_data_module.R`)
7. **Canonical data format is wide with underscore `var_<time>` naming** — Dagger's temporal auto-detection greps `_t\d+$`/`_\d+$` and would silently fail on base `reshape()`'s dot naming; both raw and pivoted frames get previews (fixes PsychoNetrix's invisible pivot). (`03_data_module.R`)
8. **`pickerInput` (PsychoNetrix's pattern) is the one variable-selection paradigm** — the only rule-3-compliant one of the three; dual-listbox and drag-and-drop dropped. (`04_varselect_module.R`)
9. **Pastel1-derived defaults replace all three source palettes** — rainbow(), sequential "Blues", and hardcoded hex arrays all dropped per rule 9. (`02_palette.R`)
10. **Dagger's unsigned arcs map onto the positive-edge picker only** — bnlearn arcs have no sign; we don't fake one. (`02_palette.R`, `05_appearance_module.R`)
11. **Figure export re-plots a stored zero-arg closure into a fresh device** — gives BootSON export it never had, adds SVG everywhere. (`05_appearance_module.R`)
12. **Package guard installs missing-only, never updates** — with the version ledger recorded as pipeline step 1 and in the script header (rule 0). (`00_packages.R`)

## Correctness wiring (network-psychometrics skill)

13. **Bootnet bootstraps are not opt-in** — "Estimate & validate" runs estimation → nonparametric boot → case-dropping boot as one recorded pipeline; centrality is locked until they complete. (`10_tab_bootnet.R`)
14. **CS gate is per-statistic and hard** — CS < 0.25 disables the ordered centrality plot; values remain visible only alphabetically under a red banner (transparency without ordering claims); 0.25–0.5 shows a caution band. (`10_tab_bootnet.R`)
15. **Gate keys on the leading interpreted measure** — expected influence when negative edges exist, strength otherwise. (`10_tab_bootnet.R`)
16. **EI force-included and listed first whenever any edge is negative** — strength misleads when opposite-sign edges cancel; both computed transparently as `colSums` on the weights matrix, not via helpers whose defaults may drift. (`10_tab_bootnet.R`)
17. **`corMethod` default flipped to `cor_auto`** — fixes BootSON's Pearson-by-default bug for Likert items; estimator/corMethod/γ pinned in every call and every fragment. (`10_tab_bootnet.R`)
18. **DAG: the bootstrap-averaged network is the reported model** — the single algorithm run is demoted to "exploratory fit", commented out in the export. (`11_tab_dagger.R`)
19. **Equivalence-class audit runs before any `cextend()`** — cextend manufactures directions; a caveat computed after it would understate what the data can't identify. (`11_tab_dagger.R`)
20. **CPDAG caveat is a named-arcs UI banner, not boilerplate** — lists exactly which arcs have unidentified direction, plus bootstrap direction-confidence (≈0.5 = undecidable) in the arc table. Dagger's "gold standard for causal discovery" copy is not carried over. (`11_tab_dagger.R`)
21. **psychonetrics estimator resolved before fitting** — the string "default" never reaches `runmodel`, the log, or the export; "auto" resolves to FIML when NAs present (skill prefers FIML), ML otherwise, and the resolved value is what's recorded. (`12_tab_psynet.R`)
22. **Ordinal+ML guard added** — warns when ordered/dichotomous data meet an ML-family estimator (the failure mode PsychoNetrix missed). (`12_tab_psynet.R`)
23. **NCT design (paired vs independent) is an explicit required radio** — `paired=` is always passed; paired path validates subject alignment (ID-matching UI is a stage-3 stub with a hard equal-n guard until then). (`12_tab_psynet.R`)
24. **NCT edge p-values corrected (Holm default)** — passed natively via `p.adjust.methods` when the installed NCT supports it, else applied manually to `einv.pvals`; significance flags computed on corrected values only. (`12_tab_psynet.R`)
25. **NCT estimation pinned to the bootnet tab's own settings reactive** — one source of truth for default/corMethod/γ, so the comparison can never silently use different settings than the networks on screen; the log also states "a null NCT is not proof of equality". (`12_tab_psynet.R`, `app.R`)

## Stage-3 boundaries

All mechanical work is marked `TODO(stage3-...)` in code: file-loader dispatch + loader UI, pickerInput/colourpicker/slider wiring, per-figure download buttons, Bridge-Symptoms + NSON port, Dagger constraint builder/split/temporal-graph port, psychonetrics non-GGM families + lambda editor port, NCT subject-ID matching UI. **Stage 3 must not alter anything above the TODO lines in the three tab files or in the recorder.**

## Post-Stage-3 user change requests (2026-07-05)

26. **Bootnet bootstraps made on-demand** (user override of decision #13) — estimation is now instant; a separate "Run bootstrap validation" button runs the nonparametric + case-dropping bootstraps. **The correctness gate is unchanged**: centrality and edge-CI panels stay locked (`validate()` refusal) until validation has run for the *current* network, and re-estimating drops the stale `bootnet_stability` fragment and nulls `rv$cs` so the gate re-locks. Only the *timing* moved, not the requirement. (`10_tab_bootnet.R`)
27. **Node size can scale by column mean** (rule-5 extension) — appearance module gains a "Scale node size by column mean" toggle + min/max range slider; both network tabs map lowest mean → smallest node, highest → largest, via a shared `scale_vsize_by_mean()` whose exact definition is emitted into the exported script. (`05_appearance_module.R`, `10_tab_bootnet.R`, `11_tab_dagger.R`)
28. **All BootSON estimators restored to the bootnet tab** — EBICglasso, ggmModSelect, pcor, cor, huge, **IsingFit** (binary, with a validity check), **mgm** (mixed), and **relimp** (relative-importance, plotted directed & non-negative). Each pins its own settings explicitly and records its own fragment. (NB: this is the *observed-variable* relative-importance network; PsychoNetrix's separate *latent* RI network (Johnson/LMG) remains a psychonetrics-tab port.) (`10_tab_bootnet.R`)
29. **DAG blacklist/whitelist constraint builder** — FROM×TO cross-product accumulation into blacklist/whitelist tables, threaded into `algorithm.args` so they apply to **every** bootstrap replicate (not just a single fit); a whitelist that forces a cycle is caught and surfaced. Constraints are baked as literal `data.frame`s into the exported script. (`11_tab_dagger.R`)
30. **DAG advanced arc metrics (direction × strength)** — arc width/labels can show boot.strength's `strength` (P present), `direction` (P orientation | present), or their `combined` product; the arc table always lists all three so nothing is hidden. (`11_tab_dagger.R`)
31. **Layout options** — bootnet tab: spring / circle / EGA-communities / PCA; DAG tab: spring / circle / hierarchical-Sugiyama. EGA is used for community colouring/layout only — the reported network is always the pinned estimator result. (`10_tab_bootnet.R`, `11_tab_dagger.R`)
