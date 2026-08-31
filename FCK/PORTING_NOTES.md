# Porting notes: how `WaPaa1_3.R` + `CIRCAREG.R` became `FCK/`

The brief was: one integrated app, all the analyses and outputs conforming to
the two existing apps, with the data (pre)processing and smoothing shared.
This is the record of exactly what that meant in practice.

## 1. What the two apps actually shared

Comparing every input and output id in the two source files:

* **25 input ids and 10 output ids collided.** All but two of the collisions
  were the *same step implemented twice* — `datafile`, `header`,
  `data_format`, `load_data`, `apply_selection`, `generate_sample`,
  `smooth_method`, `smooth_factor`, `n_basis`, `constrain_bounds`,
  `min_bound`, `max_bound`, `apply_smooth`, `run_gam_reml`,
  `run_reml_profile`, `run_cv_analysis`, `cv_k_folds`, `lambda_min_exp`,
  `lambda_max_exp`, `n_lambda`, `diag_subject_id`, `diag_subject_type`,
  `use_diagnostic_lambda`, plus `var_select_container`, `data_preview`,
  `raw_data_plot`, `data_plot`, `gam_reml_summary`, `reml_profile_plot`,
  `cv_curve_plot`, `smoothing_comparison_summary`, `diagnostics_available`.
  Those are the shared pipeline, and they now exist once.
* **Two collisions were genuinely different analyses** that happened to share
  a name: `run_pairwise` and `pairwise_correction`. WaPaa's belong to the
  post-hoc tests after functional ANOVA; CIRCAREG's belong to pairwise
  comparisons of cosinor parameters. Both tabs are kept; the cosinor one is
  prefixed `hp_`.
* **`export_scores_csv`** meant fPCA scores in WaPaa and FoSR beta
  coefficients in CIRCAREG.

The two smoothing implementations turned out to be the *same algorithm*:
WaPaa's own comments say "Applying smoothing using CIRCAREG method" and
"CIRCAREG EXACT METHOD". So sharing the smoothing step required no statistical
choice between them, only taking the union of their features.

## 2. What is hand-merged, and what is verbatim

**Hand-merged** (the shared pipeline — `FCK/ui/10_import.R`,
`ui/20_preprocess.R`, `ui/90_export.R`, `server/00_state.R`,
`server/10_import.R`, `server/20_smoothing.R`): every one carries a header
naming the source line ranges it is a union of.

**Carried verbatim by line range** — everything else. `tools/port_fck.py`
slices it out of the two source files, which stay in the repository root as
the reference. Each generated file names its provenance in its header, e.g.

```
# PORTED VERBATIM by tools/port_fck.py
#   WaPaa1_3.R lines 4425-6227  (group UIs, fANOVA (between + repeated measures))
#   WaPaa1_3.R lines 6524-6976  (post-hoc pairwise outputs)
```

`tools/port_fck.py` is a **one-time port script kept for provenance**. `FCK/`
is the source of truth from here on; re-running the script only rewrites the
files in its manifest and never touches the hand-written ones. Every edit it
makes to ported code is anchored on an exact source string, so an upstream
change turns into a hard error instead of a silent mis-port.

## 3. Coverage

Measured against the two sources:

| | outputs | helper functions | inputs read |
|---|---|---|---|
| `WaPaa1_3.R` | 91 / 91 | 26 / 26 | 95 / 96 |
| `CIRCAREG.R` | 65 / 65 (7 renamed) | 21 / 21 | 75 / 76 (3 renamed) |

The two "missing" inputs are the two variable pickers that were deliberately
merged into one step: WaPaa's `sel_group_vars` and CIRCAREG's `sel_func_vars`
(see 4.1). Nothing else from either app was dropped except one stub (4.10).

## 4. Deliberate behaviour changes

Each of these was a choice; everything not listed here behaves as it did.

**4.1 One variable-selection step.** WaPaa asked for "group variables",
CIRCAREG for "scalar variables". They are the same columns used two ways, so
there is now one picker (`sel_cov_vars`) that feeds both `values$covariates`
(original types — FoSR/SoFR/cosinor predictors and responses) and
`values$group_variables` (factor copies — fANOVA, cluster composition, cosinor
group tests). The curve columns are picked once as `sel_data_vars`.
Consequence: an Excel sheet chosen for a cosinor analysis is equally available
to the fPCA, which was not possible before.

**4.2 Separator.** CIRCAREG made the user choose; WaPaa sniffed the first line.
Both are offered, with sniffing as the default, so files that worked in either
app still load.

**4.3 Excel sheet selection** (CIRCAREG only) now applies to every analysis.

**4.4 Smoothing** is WaPaa's implementation — which captures per-subject EDF
and GCV, reports relative RMSE as a percentage of the data range, allows a
different basis count for automatic and manual modes, and builds the `fd`
object on a 0–1 range because fPCA, warping, fANOVA and clustering all assume
that range — **plus CIRCAREG's cyclic option**: ticking "Is data cyclic?"
switches to a Fourier basis with `min(n_time, 13)` functions and lowers the
minimum-points threshold from 4 to 3, exactly as CIRCAREG did. When cyclic, the
0–1 representation is built with a Fourier basis too, so the periodicity
survives into the fPCA.

**4.5 Smoothing diagnostics** are WaPaa's, which is a strict superset of
CIRCAREG's: it adds group-stratified CV folds and the GCV-vs-n-basis sweep, and
is otherwise line-for-line the same code. *Known divergence:* inside the
cross-validation WaPaa fixes the basis at `min(20, n_time - 2)` where CIRCAREG
used the user's `n_basis`. WaPaa's is kept so its numbers do not move; if you
smooth with a very different basis count, read the CV curve with that in mind.

**4.6 Sample data.** The two generators produced different datasets (WaPaa: 100
points on 0–1 with group structure; CIRCAREG: 24 hourly points with covariates
and a binary outcome). One dataset has to serve every tab, so the merged
generator produces a 24-hour circadian set: 50 subjects, hourly columns
`00:00`–`23:00`, group structure with real MESOR/amplitude/acrophase
differences, `Age`/`Sex`/`Score` covariates, and a binary `Outcome` driven by
the curve level. It exercises all seven analyses.

**4.7 Re-smoothing clears downstream results.** *(new)* In the source apps a
stale fPCA or cosinor fit could survive a re-smooth. With one smoothing step
feeding seven analyses that would be much easier to miss, so applying smoothing
(or a new file, or a new variable selection) clears the analysis results and
they have to be re-run.

**4.8 Harmonic time values.** The cosinor tab detects clock times from the
column names itself; the shared import now does the same thing once, into
`values$time_numeric`. Rather than have two detectors disagree, the cosinor
time selector gained one extra option — "Use shared times detected at import".
It is additive: the default is still `_index_`, so the tab behaves exactly as
before unless you pick the new option.

**4.9 Landmark plot NULL guard.** *(bug fix)* WaPaa's landmark plot branches on
`input$landmark_target` and `input$selected_subject`, neither of which its UI
ever creates. Both are `NULL`, so `NULL == "mean"` is `logical(0)` and `||`
raises "invalid length zero argument" on R ≥ 4.3 — the plot errored instead of
drawing the mean curve. The test is reordered so the `NULL` check comes first;
behaviour is unchanged if those inputs are ever supplied.

**4.10 One R-code export.** CIRCAREG's "Download Reproduction R Code" button
wrote a single placeholder comment (`"# R Code Generation logic here (omitted
for brevity in single file app)"`). It is dropped rather than shipped as a
working-looking button; WaPaa's real generator is the app's code export.
CIRCAREG's FoSR coefficient export survives as `export_fosr_coefs_csv`.

**4.11 Optional packages no longer block startup.** Both source apps called
`install.packages()` at load and stopped if any failed. With seven analyses
behind one launcher, a failed `reticulate` install would have taken the cosinor
tabs down with it. The twelve packages the app cannot run without are still
required; the seven that back a single feature each are loaded if present and
reported if not.

**4.12 Tab id.** CIRCAREG's pairwise tab used `tabName = "pairwise"`, which
WaPaa's post-hoc tab already owns; the cosinor one is now `harm_pairwise`.

## 5. Rename table

| source | source app | merged app |
|---|---|---|
| `sel_group_vars`, `sel_cov_vars` | W, C | `sel_cov_vars` (one picker) |
| `sel_data_vars`, `sel_func_vars` | W, C | `sel_data_vars` (one picker) |
| `values$gam_reml_result` | C | `values$gam_reml_fit` |
| `values$reml_profile_result` | C | `values$reml_profile` |
| `values$cv_result` | C | `values$cv_results` |
| `values$pairwise_results` | C | `values$hp_pairwise_results` |
| `values$pairwise_param` / `_correction` | C | `values$hp_pairwise_param` / `_correction` |
| `run_pairwise`, `pairwise_correction`, `pairwise_param`, `pairwise_show_ci`, `pairwise_show_effect_size` | C | `hp_run`, `hp_correction`, `hp_param`, `hp_show_ci`, `hp_show_effect_size` |
| `pairwise_results`, `pairwise_plot`, `pairwise_matrix`, `pairwise_matrix_help` | C | `hp_results`, `hp_plot`, `hp_matrix`, `hp_matrix_help` |
| `export_pairwise_results`, `export_pairwise_plot` | C | `hp_export_results`, `hp_export_plot` |
| `export_scores_csv` | C | `export_fosr_coefs_csv` |
| `tabName = "pairwise"` | C | `tabName = "harm_pairwise"` |

The first three `values$` renames are smoothing diagnostics: the merged app
runs WaPaa's version of that section, so CIRCAREG's names disappear with its
duplicate code.

## 6. Known gaps

* **The R-code export covers the WaPaa pipeline only.** `export_code` writes a
  script that reproduces the import, smoothing, fPCA, fANOVA and clustering
  steps. The cosinor, FoSR and SoFR fits are *not* written into it — CIRCAREG
  never had a generator to port. Their parameters export as CSV instead.
* **Two time-detection paths still exist.** The cosinor tab keeps its own
  column-name parser alongside the shared one (4.8). They agree on the common
  formats; the shared one is one click away if they ever disagree.
* **`group_summary` and `group_preview_plot`** are defined in the server but no
  UI places them. That is true in WaPaa as well — dead code carried across
  rather than silently deleted.
* **The CV basis-count divergence** in 4.5.
* **Not run end-to-end here.** `tests/smoke_test.R` verifies that every file
  parses, the whole UI renders, every tab is reachable, no output is defined
  twice and all 16 server files register in one environment. It cannot verify
  the statistics: `fda`, `refund`, `rmfanova`, `fda.usc` and `shinyWidgets`
  were not installable in the environment this port was done in (CRAN was
  blocked by egress policy), so the analyses themselves have not been executed
  since the merge. Run each tab once against a known dataset before trusting
  the numbers.
