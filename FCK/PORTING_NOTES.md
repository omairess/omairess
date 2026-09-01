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
is otherwise line-for-line the same code. One line is changed: inside the
cross-validation WaPaa fixed the basis at `min(20, n_time - 2)`, where CIRCAREG
used the user's `n_basis`. CIRCAREG's behaviour wins here — the CV curve exists
to recommend a smoothing factor for the smoothing you are about to run, and a
basis count unrelated to yours makes the recommendation mis-targeted.

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

**4.8 Harmonic time values, and the two meanings of "time".** WaPaa has two
notions of time that are easy to confuse:

* `extract_time_values()` returns `1:n_time`. It does **not** read the column
  names — its own comment says "Use sequential numbering to maintain
  chronology". This is what fills `values$time_numeric` and what every ported
  plot uses as its x coordinate.
* `extract_hour_from_colname()` is the real parser (`Base9h` → 9,
  `Base7h30` → 7.5, `KSS_9u_dag1` → 9, `X8.25` → 8.25) and is used only to
  build axis *labels*.

Both are left exactly as they are, so every ported plot keeps its original
coordinates. `FCK/server/03_helpers_clock.R` adds a third, explicitly named
thing: `values$time_clock`, real clock hours reported **only when they can be
trusted** — every column must parse to a finite hour in `[0, 24)`, which
rejects the `T1…T30` case where the parser's "any number" fallback would
otherwise hand back column indices dressed up as hours.

On that basis the cosinor time selector gains one extra option, "Use shared
clock times parsed at import". It is additive: the default is still `_index_`,
so the tab behaves exactly as before unless you pick it, and it refuses with a
message rather than guessing when no clock times could be parsed.

*(This corrects an earlier version of the merge in which that option read
`values$time_numeric` — i.e. it fed `1, 2, 3, …` to the cosinor fit while its
label promised clock times. Harmless for evenly-spaced hourly columns, wrong
for the unequally-spaced designs the option exists to serve.)*

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

**4.13 Smoothing can use real clock times.** *(new, opt-in)* Both source apps
smoothed against the column index (`t <- 1:n_time`), which is correct only when
the columns are evenly spaced in real time. CIRCAREG's own help text
acknowledges the problem — it tells you to enter times manually "if your
measurements are unequally spaced (e.g. hourly during day, 2-hourly at night)"
— but that escape hatch existed for the *cosinor fit only*. The smoothing
underneath it, and therefore the fPCA, fANOVA and clustering built on the
smoothed curves, still treated a 2-hour gap as one step exactly like a 1-hour
gap.

A checkbox on the smoothing tab, **"Space time points by their real clock
times"**, switches the smoothing argument to elapsed hours
(`fck_cumulative_hours()`, unwrapped across midnight: 20, 22, 2, 8 → 0, 2, 6,
12). The basis range follows, a cyclic Fourier basis is then given
`period = 24` explicitly rather than fda's default of "however long the
recording happened to be", and the 0–1 representation the fPCA family reads is
rescaled *preserving the spacing* — otherwise the `fd` object would put the
columns back on an even grid and undo the whole thing.

It is **off by default**, so existing results reproduce exactly. When on, the
roughness penalty becomes per hour rather than per column, so the smoothing
factor needs re-tuning; the app says so when you apply it, and the fit summary
now reports which time axis produced the numbers. If the column names yield no
usable clock times the app warns and falls back to the column index rather than
silently doing something else.

**4.14 "No smoothing" now means no smoothing.** *(bug fix)* The fPCA / warping
/ fANOVA / clustering family needs an `fd` object, and building one always
means projecting onto a basis. The source apps wrote that projection in four
places with two different sizes, and in two of them it fired without being
asked:

* choosing **"Raw data (no smoothing)"** still built `fd_obj` on
  `min(20, n_time - 2)` basis functions — on 24 hourly columns that is a real
  smooth, applied under a control labelled *no smoothing*;
* running fPCA or fANOVA **without visiting the smoothing tab** did the same
  thing, silently, from inside the analysis tab.

There is now one rule, in `FCK/server/04_helpers_fd.R`: no smoothing means an
*interpolating* basis (`nbasis = n_time`), so `smooth.basis` is a square system
that reproduces the observed values at the observed times. The analysis tabs
call `fck_ensure_fd_obj()`, which builds that representation **and says so** in
a notification rather than choosing a smoothing basis on the user's behalf. The
penalised smoothing on the preprocessing tab is untouched and remains the only
place a roughness penalty is applied.

To be explicit about what this does *not* claim: an `fd` object cannot avoid
committing to some behaviour between the measurement points. What it can avoid
is changing the values *at* them, and now it does.

**4.15 The exported R script covers every family.** *(C1)* WaPaa's generator
covered its own analyses; CIRCAREG had none. Three sections are added —
harmonic/cosinor, FoSR, SoFR — driven by settings each fit records about itself
(`fck_settings`) rather than by reading the widgets at export time, which would
emit whatever they happen to say now instead of what produced the stored model.

The cosinor section emits the app's own `fit_cosinor()` through `deparse()`
rather than a hand-written re-implementation. A re-implementation is a second
copy of 300 lines of fitting logic that can drift from the one that produced
the numbers; dumping the real function cannot. `tests/codegen_test.R` drives
the generator with a stub result for every family and checks the output is
valid R.

**4.16 Sessions can be saved and reopened.** *(C2)* Neither source app could
be: close the browser and the import, the smoothing and every fitted model were
gone. `FCK/server/91_session.R` writes one `.rds` holding the whole `values`
bus, the analysis-defining settings, and the exact package versions. On
restore it reports what came back and **warns if any package version has moved
since**, because an `fda` or `mgcv` update can shift the numbers a restored
model is sitting next to. Cosmetic display controls come back at their
defaults; only the settings in `RESTORE_INPUTS` are pushed back into widgets,
because restoring every input generically means guessing each widget's message
format and failing silently on the ones it gets wrong.

**4.17 The group-structure views are visible.** WaPaa defined
`group_summary` and `group_preview_plot` but no tab ever placed them — dead
code in the source app. They are now a collapsible box on the Data Import tab,
next to the selection that produces them.

**4.18 Filled points are visible, and extrapolation is off by default.**
*(new)* Smoothing fills every gap — a subject measured at 8 of 20 times comes
out of the smoother with 20 values. Both source apps did this correctly and
neither showed which values were invented.

The distinction that matters is not observed-vs-missing but *where* the missing
value sits. Filling **between** two of a subject's own measurements is what
smoothing is for. Filling **beyond** their first or last measurement is not:
`eval.fd()` is evaluated across the whole grid, and a B-spline carried past its
data follows its end polynomial rather than levelling off.

On a real staggered sleep-deprivation protocol (13 participants x 7 cycles, 20
two-hourly columns spanning 06:00 day 1 to 20:00 day 2) the split is **46 %
observed, 22 % interpolated within range, 32 % extrapolated beyond it** — a
median of 12 h and a maximum of 32 h invented per row, with 7 of 85 rows more
than half filled in. One row has 3 measurements and 16 extrapolated points.
That is not an edge case; it is what a staggered protocol looks like.

So:

* `values$fill_status` records, per cell, whether it was observed, interpolated
  or extrapolated (`server/05_helpers_missing.R`), computed from the raw data
  before anything is filled.
* A **"Missing data & filled points"** panel on the smoothing tab shows the map
  (subjects x time, three states, sortable), a per-subject table sorted worst
  first, and a single-subject inspector plotting the curve with its measured
  points solid, interpolated points hollow, and extrapolated points crossed.
  It downloads as a CSV alongside the smoothed curves.
* **Extrapolation is off by default**: values beyond a subject's observed range
  are held flat at their first/last measurement rather than following the
  fitted polynomial. This is a deliberate change from both source apps, which
  always extrapolated. Flat-held values are still invented and the map still
  marks them; the checkbox restores the original behaviour.
* An optional **minimum measured points per row** at import drops rows that are
  mostly reconstruction. It defaults to 0 (keep everything) and is applied where
  the frame is rebuilt from `raw_df`, so lowering it and pressing Confirm again
  brings the rows back.

`tests/missing_status_test.R` pins the classification, including the staggered
pattern above; mislabelling an extrapolated point as interpolated would be
invisible, since the map would simply reassure you about a value nothing
constrains.

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

* **The exported script is checked for syntax, not for equivalence.**
  `tests/codegen_test.R` proves it is valid R and that every family
  contributed a section; nobody has yet run an exported script end to end and
  compared its numbers against the app's. The cosinor section emits the app's
  real fitting function, so it should agree exactly; the FoSR GAM section
  reproduces the model but derives its coefficient curves the same
  approximate way the app does (prediction contrasts), which is exact for
  numeric predictors and approximate for factors.
* **The FoSR bootstrap has no seed in the app.** Its bands differ slightly run
  to run. The exported script fixes `set.seed(1)` so the script is
  reproducible, which means the script's bands will not match any particular
  run of the app exactly.
* **Two time-detection paths still exist.** The cosinor tab keeps its own
  column-name parser alongside the shared one (4.8). They agree on the common
  formats; the shared one is one click away if they ever disagree.
* **Real-time spacing is opt-in, not automatic.** The Data Import tab now says
  when your columns are unevenly spaced, but it will not switch the smoothing
  for you — that would change results for anyone who re-ran an old analysis.
* **Real-time smoothing has not been run end-to-end.** Its parsing and midnight
  unwrapping are covered by `tests/clock_helpers_test.R`, but the fda calls
  themselves have not been executed (see below).
* **`group_summary` and `group_preview_plot`** are defined in the server but no
  UI places them. That is true in WaPaa as well — dead code carried across
  rather than silently deleted.
* **The CV basis-count divergence** in 4.5.
* **Not run end-to-end here.** `tests/smoke_test.R` and
  `tests/clock_helpers_test.R` verify that every file
  parses, the whole UI renders, every tab is reachable, no output is defined
  twice and all 16 server files register in one environment. It cannot verify
  the statistics: `fda`, `refund`, `rmfanova`, `fda.usc` and `shinyWidgets`
  were not installable in the environment this port was done in (CRAN was
  blocked by egress policy), so the analyses themselves have not been executed
  since the merge. Run each tab once against a known dataset before trusting
  the numbers.
