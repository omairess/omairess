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

**4.19 The post-hoc tests compared the wrong variable.** *(bug fix)* The omnibus
functional ANOVA runs on `input$fanova_group_var` (via
`get_fanova_group_labels()`), optionally restricted to a subset of its levels on
an fd object subset to match. The post-hoc tests then called
`perform_pairwise_comparisons()` with `values$group_labels` — the **primary**
grouping variable, i.e. whichever scalar variable happened to be selected first
at import.

With two or more scalar variables selected, the omnibus tested one variable and
its "post-hoc" tests silently tested another, on the full unfiltered curve set.
Nothing in the output said so, which is why `tests/posthoc_source_test.R` pins
it rather than leaving it to inspection.

Resolution now happens in one place (`FCK/server/06_helpers_posthoc.R`), driven
by a control on the post-hoc tab:

* **Follow the functional ANOVA** (default) — the labels, the levels, the
  curves and the design the omnibus actually used, read back from
  `values$fanova_results` (which now also carries `fd_used` and `group_var`).
  They cannot drift apart, because there is only one source.
* **Choose a variable and design here** — including a **within-subjects
  (paired)** comparison, which previously required the omnibus itself to have
  been a repeated-measures ANOVA. A comparison on a variable the omnibus was
  never run on is a *new family of tests*, not post-hoc ones; the app labels it
  that way rather than letting it borrow the omnibus's authority, and the
  correction is stated to cover that family only.

The resolver refuses rather than guesses: a single-level variable, labels that
do not line up with the curve count, or a paired design where no subject appears
at more than one level all stop with a message naming the problem.

**4.20 Repeated-measures columns were read from unfiltered rows.** *(bug fix)*
The repeated-measures path read its subject-ID and factor columns from
`values$uploaded_data` — the raw imported frame, **before** the import step
drops all-missing rows (and, in this app, rows below the minimum-measured-points
threshold). Any dropped row shifted every pairing, and nothing errored: the
lengths only had to be plausible for the test to return numbers.
`fck_rm_column()` prefers the row-aligned `values$covariates` and, when it must
fall back to the raw frame, checks the length and stops with an explanation
instead of mispairing.

**4.21 A circular density of the acrophases.** *(new)* The existing polar plot
draws one marker per subject at (amplitude, acrophase), which overplots badly at
85 rows and shows where individuals are rather than where the sample sits. Tab
**2b. Acrophase Density** puts a von Mises kernel density on the same circle —
the circular analogue of a KDE, so mass near midnight wraps instead of being
split between the two ends of a histogram.

The circle is oriented as a clock face: **noon at the top, midnight at the
bottom**, hours running clockwise, so the upper half is daytime (06:00–18:00)
and the lower half night, with the night sector shaded (the dusk/dawn bounds are
adjustable). In plotly's angular convention that is
`theta(h) = (270 - 15h) mod 360`.

The radius carries the **fitted harmonic-regression curves** (the default — the
curves from tab 1 wrapped onto the clock, built from the same `mean_coefs` and
the same band, so it is that plot in polar coordinates rather than a second
opinion about the same data), the acrophase density, or the **signal itself
averaged over the clock** (the smoothed curves wrapped onto one 24 h face). A recording
longer than the period visits the same clock time twice — a 38 h protocol hits
06:00 on both days — so the profile mode either averages the days at each clock
time or draws one ring per day; it never silently folds them together. The shape
is filled, and an adjustable **inner radius** keeps a low stretch of the clock
from collapsing to a point, which is what turns a ring into petals; the note
says plainly that a non-zero inner radius breaks proportionality between radius
and value, and 0 restores it.

Bandwidth defaults to Taylor's (2008) circular plug-in rule, **capped at
period/12**, and can be set by hand in hours. The cap is not cosmetic: the rule
estimates concentration from R̄, a *global* quantity that collapses for
multimodal data — symmetric bimodal acrophases have R̄ ≈ 0 however tight each
mode is — so the uncapped rule flattened the ring into a featureless disc
exactly where there was structure to see (measured max/min 1.1 on tight bimodal
data). Deriving the concentration from higher trigonometric moments was tried
and rejected: it rescued the multimodal cases but grew lobes in genuinely
uniform data (max/min 6.8), and a density plot that manufactures structure is
worse than one that misses a subtle mode. The cap behaves across every shape
tested (tight bimodal 14.8, uniform 1.7) and cannot manufacture anything. The
note reports the max/min contrast and, when the ring is nearly round, says that
this is the data rather than a broken plot.

Fixed while measuring this: `besselI` overflows to `Inf` for large κ, so the
rule fell through to its fallback — the *widest* bandwidth — for the *most
concentrated* samples, which is backwards. The exponentially scaled forms cancel
exactly (`I2(2κ)` and `I1(κ)²` both carry `exp(2κ)`), so very tight acrophases
now get 0.25 h rather than the ceiling. Amplitude weighting is optional and off by default: unweighted, a
barely-detectable rhythm counts as much as a strong one even though its
acrophase is mostly noise. The mean direction is drawn with length = R̄, and the
note under the plot states the bandwidth, n, the circular mean, R̄ and what that
concentration means.

Colour follows the reference categorical palette in fixed slot order. Overlapping
density rings are an all-pairs comparison, where only the first three slots clear
the colour-vision gates (worst pair ΔE 9.2 deutan, 24.0 normal-vision), so line
style distinguishes series as well as hue and the note says so past three groups.
Non-24 h periods are plotted as phase within the period, with the day/night
framing explicitly disclaimed.

`tests/circular_density_test.R` pins the orientation (a sign slip mirrors the
clock, and a mirrored clock still looks like a clock), the wrapping, that the
density integrates to 1, that a very small bandwidth does not overflow, and the
weighting behaviour.

Two things about the fit mode are worth stating, because both look like bugs
and are not:

* **A cosinor in polar coordinates is a limaçon.** `r = MESOR + amplitude *
  cos(angle - acrophase)` traces an *off-centre ring*, widest toward the
  acrophase — not a lobed blob. The offset is the rhythm. The note says so.
* **A trend is not drawn.** A linear, log or saturating trend is not periodic,
  so 08:00 on two different days would sit at the same angle with different
  values and the ring would not close. The polar plot shows the rhythm and says
  where the trend is.

The second needed care in the code. `predict_from_coefs()` lays coefficients out
as `c(mesor, [trend...], beta_cos_1, beta_sin_1, ...)` and locates the harmonics
at an offset that depends on `trend_type`, so calling it with
`trend_type = "none"` to drop a trend would read the **trend** coefficient as the
first harmonic — a silent, entirely plausible-looking wrong curve.
`fck_rhythm_from_coefs()` keeps the real `trend_type` for indexing and simply
never adds the trend term. `tests/circular_density_test.R` pins that it equals
the fitted curve minus its trend, and explicitly that it is *not* what dropping
`trend_type` would give.

Also pinned there: a known amplitude and acrophase survive the round trip
(a curve built to peak at 03:30 peaks at 03:30 on the clock), and a missing
harmonic coefficient is skipped rather than poisoning the whole curve.

**Night is drawn as a gradient**, fading in at dusk, deepest halfway through the
window and fading out at dawn, over a default window of 23:00–07:00. A hard grey
half-circle split at 06:00/18:00 asserts an edge that night does not have; a ramp
carries the same information without the assertion. It is drawn as an *annulus*
rather than wedges from the centre — wedges converge to a point, so the colour is
most intense at the origin, the one place on the plot that carries no time at
all, and it tints the data fills worst where they overlap. The flat block is
still available, as is none.

**Every trace of a series carries `legendgroup`.** Splitting fill from line is
what lets the lines sit above every fill, but without a legend group a legend
click hides only the trace it names: the line vanishes and its fill stays behind,
which reads as "that curve went faint" rather than "that group is off".
`tests/circular_density_test.R` checks this at source level, since verifying it
in the render would need a browser.

The radial axis is **labelled in the response's own units**, with the tick radii
computed through the same mapping the data went through, so a label at 60 sits
exactly where a fitted value of 60 is drawn. Angular ticks default to every
hour. Uncertainty is drawn as **dotted edges** by default rather than a shaded
ribbon: a ribbon hides the ring wherever two groups overlap, and dotted edges
stay legible over a fill and over each other (shaded, both, or none are also
available). Fills are drawn for every ring first and the lines on top of all of
them — ring-by-ring, the second group's fill lies over the first group's line
and hides exactly the comparison the plot exists to make.

A radial-baseline control chooses between measuring the radius from the smallest
value plotted (stretches the shape, small differences visible) and from zero
(proportional, but flattens a rhythm whose MESOR is far from zero). The note
states which is in force, since neither is right for every reading.

**4.22 The preprocessing plot drew a fixed 50 curves.** *(new control)* WaPaa
capped its curve plot at `min(n_subj, 50)` and said nothing about the rest, so on
85 rows a third of the data was silently absent from the picture people use to
judge their smoothing. The count is now a control (10 / 25 / 50 / 100 / 250 /
All), defaulting to 50 so nothing moves for an existing analysis.

It had to change in **two** places. The click-to-select handler recomputes the
same cap to work out which trace was clicked, so leaving that at 50 would mean
clicking the 60th curve selects nothing or the wrong subject.

**4.23 The polar tab disagreed with the two plots it mirrors.** *(bug fix)*
Reported from screenshots: the same clock time read 86.6 on the Harmonic
Regression Fit, 47.5 on the polar fit ring and 21.4 on the averaged-signal
ring. Three numbers, one model. Two separate causes, both mine, both silent —
each produced a plausible circle rather than an error.

*The fit ring dropped the trend.* `fck_rhythm_from_coefs()` computed the
coefficient offset for the trend parameters but never added the trend TERM,
so it drew MESOR + harmonics while `predict_from_coefs()` on tab 1 drew
MESOR + trend + harmonics. The two therefore differed by exactly the trend at
that time — about 39 units at 03:00 of the second day on the reporter's data,
which is most of the signal when the trend dominates. The omission was
deliberate and documented (a trend is not periodic, so 08:00 on two different
days is one angle with two values and the ring cannot close), but the plot
never said so and labelled its radial axis "fitted response". A design
decision the reader cannot see is indistinguishable from a bug.

The fix keeps both curves and makes the tab say which one is on screen.
`fck_rhythm_from_coefs()` gains `include_trend` and `t_offset` and now mirrors
`predict_from_coefs()` term for term. A new control, **Curve shown**, offers:

* *The fit over the recording* (default) — evaluated on the model's own
  absolute axis (`mod$time_vec`: 8, 9, … 30) and placed on the dial at
  `t %% 24`, so it matches tab 1 value-for-value. It is an **open arc**: it
  covers only the hours recorded, and the gap between its two ends is the
  trend. Filled back along the inner radius rather than with `toself`, so the
  uncovered hours read as empty instead of as a flat stretch of curve.
* *The rhythm only* — the periodic part, which is the only thing that closes
  into a ring. The radial axis relabels itself to `fitted response, trend
  removed` and the note states that the values are not tab 1's.

With `trend_type = "none"` the two coincide and the control is inert.

*The averaged-signal ring was rotated by the start hour.* `fck_profile_rings()`
took `fck_cumulative_hours()` — which returns **elapsed** hours from the first
column, starting at 0 — and passed it to `fck_wrap_to_clock()` as if it were a
clock time. Column 1 landed at 00:00 whatever time the recording began. On an
08:00 start the whole ring was rotated 8 h backwards: the 04:00 peak was drawn
at 20:00, and hovering 03:00 returned the value from 11:00. That is a real
result pointing at the wrong time of day, which is worse than a missing one.
It now reads the clock hour that `fck_clock_hours()` already parsed per column
and uses `cum` only for the day index.

`fck_open_path()` and `fck_band_path()` are the unsorted counterparts of
`fck_close_ring()` and `fck_band_ring()`. Sorting is right for a periodic ring
and wrong for an arc running 08:00 → 06:00 through midnight, where it would
reorder the path to 00:00 → 23:00 and draw the curve backwards through itself.

`tests/polar_agreement_test.R` covers both: the polar evaluator is checked
term-for-term against a transcription of `predict_from_coefs()` under linear
and saturating trends and two harmonics; the trend-free curve is checked to
close and the fit to *not* close; and the elapsed-vs-clock distinction is
asserted directly, including that the offset between them is exactly the start
hour.

**4.24 Harmonic-regression audit.** *(bug fixes + statistical changes)* A full
audit of the cosinor pipeline — fitters, aggregation, reporting — against a
supplied brief. Seven hard bugs and six statistical deficiencies were named in
the brief; five more bugs were found in the code that were not. The complete
changelog is at the top of `server/72_harmonic.R`; the shared arithmetic moved
to the new `server/08_helpers_cosinor.R`, which carries its own.

Highlights, in roughly the order they cost the most:

* **The pooled fitted equation dropped the homeostatic term** for *every* trend
  type, not just `exp_sat`. The pooled builder read `pop$indiv_means$A_sat` and
  friends; `indiv_means` was never given a single trend parameter, so the
  branch was dead. The group builder read a different structure, which is why
  only the pooled line was wrong. One renderer now, `fck_format_equation()`.
* **The Rayleigh test ran on the amplitude-weighted resultant** (0.824, Z = 886)
  where it is defined on unit vectors (0.789, Z = 812). Both resultants are now
  returned under names that cannot be swapped.
* **The zero-amplitude F test credited the trend to the rhythm** — not in the
  brief. The numerator was the whole model's SS over the harmonics' df alone.
  This is a bigger contributor to the "95.3% significant" rate than the FDA
  smoothing is.
* **Convergence was never checked.** `nls(warnOnly = TRUE)` returns a fit at the
  iteration limit rather than erroring, so every non-converged optimisation was
  counted as a success — which is how `Successfully fitted: 1305 / 1305` sat
  next to an R² range starting at 0.060.
* **"LOOCV RMSE" in the nonlinear path was in-sample residual RMSE** — not in
  the brief — while the linear path under the same label did real LOOCV.
* **The MESOR was not a MESOR.** It is the fitted constant, now called the
  intercept everywhere including the plots, the parameter table, the CSV export
  and the pairwise comparisons. The rhythm-adjusted mean is computed by
  integration and reported under the correct name. Note the brief's own sanity
  figure was slightly off: the fitter anchors the trend at the *first
  observation*, not midnight, so the mean S-factor over [8, 30] is 0.4975 and
  the rhythm-adjusted mean is about 43.8, not the 50 an unshifted `t` would give.
* **`acrophase_time_h` omitted the `/h` divisor** — not in the brief — putting
  H2 on a 0–24 scale while the group summaries used 0–12, so the H2 group tests
  ran on the wrong scale.

Estimators deliberately changed: the Rayleigh Z, the variance decomposition
(commonality, replacing overlapping marginals), the zero-amplitude F test, the
nonlinear amplitude/acrophase SEs, and the nonlinear marginal R²_S.

`tests/testthat/test-cosinor-audit.R` implements the brief's regression values;
`tests/audit_test.R` runs that same file under a base-R shim, because testthat
is not installable here and an unexecuted test is a claim rather than a check
(298 assertions pass). `tests/report_harness.R` generates `report_OLD.txt` and
`report_NEW.txt` from synthetic data with known ground truth — the real
1305-subject file is not in this repository, so old and new are scored against
truth rather than diffed against a re-run that could not be produced honestly.

**4.25 Per-component ANOVA on fPCA / warped-fPCA scores.** *(new feature)* The
fPCA tab produced an n x k score matrix and stopped; there was no way to ask
whether groups differ on a component. `server/09_helpers_pcanova.R` (arithmetic,
testable without a session) and `server/42_fpca_anova.R` (controls, report,
plot, CSV) add it, with the section appended to `ui/41_results.R`.

What it does beyond a stack of `aov()` calls:

* **Two multiplicity families, corrected separately.** Pairwise-within-a-
  component and omnibus-across-components are different families. Correcting
  the first and forgetting the second is the usual error and the more serious
  one: k omnibus tests are what manufacture a "significant PC" by chance. Both
  have their own control.
* **The post-hoc gate runs on the across-component-adjusted p**, so a component
  that only survives before that correction gets no pairwise table lending it
  credibility. `tests/testthat/test-pc-anova.R` pins this with an effect sized
  to pass raw (p = 0.028) and fail Bonferroni over 8 components (p = 0.22).
* **Corrections offered:** none, Bonferroni, Holm, Hochberg, Hommel, BH, BY,
  Tukey HSD, Games-Howell. Tukey and Games-Howell are family-wise by
  construction, so `p.adjust` is *not* applied to them a second time — a test
  asserts that.
* **Welch by default when the variances differ**, chosen per component by
  Brown-Forsythe. With n = 654 against n = 59 the equal-variance F is not a safe
  default. Kruskal-Wallis is reported alongside for the severely non-normal case.
* **Eigenvalue separation** is shown next to every component: when two
  eigenvalues are nearly tied the split between their eigenfunctions is an
  arbitrary rotation, and a "difference on PC2" belongs to the pair.
* **Repeated curves are detected.** Scores are one row per *curve*;
  `values$subject_ids` is now captured at import and carried through every row
  filter, so the report can say that a between-groups test on repeated rows is
  anticonservative instead of silently being one.

One bug found in this module while testing it, worth recording because it is a
general R trap: `l[[j]] <- NULL` **deletes** element j and shifts everything
after it down. Building the per-component lists that way meant a single
untestable component renumbered every later one, so PC4's result would be
reported as PC3's. Fixed with `l[j] <- list(NULL)`; a test pins it.

**4.26 The environment opened up.** Sections 4.24 and earlier record that CRAN
was blocked and that the statistics had never been executed. That is no longer
true and the claims are corrected here rather than left standing:

* `minpack.lm`, `testthat`, `mgcv` and `readxl` install from the Ubuntu archive
  (`apt-get install r-cran-*`), which is reachable where CRAN is not.
* `fda` is not packaged for Ubuntu; it installs from the read-only CRAN mirror
  on GitHub (`github.com/cran/fda`) **with the `fds` dependency removed** —
  `fds` is a data package that fda's R code never calls, only its NAMESPACE
  imports. Basis construction, `smooth.basis`, `fdPar`, `eval.fd` and `pca.fd`
  all verified working. Anything touching fda's bundled example datasets will
  not.
* The full `testthat` suite now runs for real: **405 assertions pass**, and the
  base-R shim in `tests/audit_test.R` was confirmed faithful (it reported the
  same 298 for the cosinor file).
* `tests/real_data_run.R` runs the audited pipeline on the real Circaflex file
  by extracting the app's own fitters from `server/72_harmonic.R` via the parse
  tree — the shipped code, not a copy. `tests/real_fpca_anova.R` does the same
  for the component ANOVA.

**4.27 Two user-reported defects in the harmonic tab.** *(bug fixes)*

**`Error in t.test.formula: grouping factor must have exactly 2 levels`.**
`n_groups` was computed once on the full `params` frame and then reused to
choose t-test vs ANOVA inside blocks that had since taken a *subset* of it —
`params_trend <- params[!is.na(params$A_sat), ]`, and similarly for `tau`. If
that filter leaves only one group with a usable trend parameter, the branch and
the data disagree: `n_groups` says 2, the subset does not, and `t.test` stops.

It became reachable when the convergence gate (4.24 / audit 2.3) started
excluding non-converged fits. On the Circaflex file that removes 938 of 1305
subjects, and a sparse age band can lose every member that had a usable
`A_sat`. The latent bug is older than that change; that change is what exposed
it, which is mine to own.

Worth recording: the obvious diagnosis is wrong. This is **not** about a factor
keeping unused levels — `t.test.formula()` calls `factor()` on the grouping
column, which drops them itself, so an undropped 4-level factor holding 2
values works fine (verified on R 4.3.3). A test asserts that, so nobody
re-derives it. A second genuine route exists and is also pinned:
`length(unique())` counts `NA` as a group, so one real group plus missing
labels reads as two.

The fix is structural rather than another guard. Every comparison in the block
now goes through one `.group_test()` entry point that reads the number of
groups off the data it is *actually handed*, after dropping unused levels and
rows with missing values, so the branch and the test can never disagree again
whatever filter ran in between. It declines politely (with a reason) at fewer
than two usable groups or fewer than two observations in any group, instead of
stopping the whole report. The block also returns early, naming the convergence
counts, when fewer than two groups survive at all. Six hand-rolled
t-test/ANOVA pairs collapsed into it; the two-group-only Cohen's d they carried
is replaced by `fck_group_linear_test()`, which handles any number of groups.

**The fit plot showed linear model time on its axis and hover.** Pointing at
the curve reported `27.01508` instead of `03:01`. The model is fitted on
unwrapped linear time and has to be — the harmonics would not care, since
`cos(2πht/T)` is identical at t = 3 and t = 27, but a **non-periodic trend
needs 08:00 on two different days to be two different values of t**, or the
days collapse onto one point and the homeostatic rise cannot be estimated. So
linear time stays the computational axis and clock time is purely a display
transform: `fck_clock_label()` and `fck_clock_ticks()` in
`server/08_helpers_cosinor.R`, applied to the axis ticks and to the hover text
of every trace. Nothing in that path feeds a fit.

Four things were wrong with the old labelling, all fixed: labels appeared only
when the recording happened to wrap past the period, so a within-day recording
got bare numbers; the label was `paste0(t %% period, ":00")`, which renders a
half-past tick as `8.5:00`; there was no zero padding and no day marker, so
08:00 on the first day and on the second were indistinguishable; and the hover
was never converted at all. Midnight crossings now get a dotted vertical rule,
and the y axis picks up the DV name and units when they have been declared.

**4.28 A recording longer than 24 h on a 24 h dial.** *(display fix + new
option)* Reported as "the harmonic fit and the circular plot of the harmonic
fit don't match again".

**They do match.** `tests/testthat/test-polar-laps.R` evaluates both paths from
the same coefficients and asserts equality at every point, per group: the
largest difference is 0 to machine precision, and the peak and trough clock
times agree. What differs is what a recording longer than one period *looks
like* on a 24 h dial — it **laps**, passing over each clock hour more than
once at different values, so the drawn shape resembles nothing in the 2-D plot
even though the numbers are identical. Three changes follow from that:

* The note now prints the agreement as a **computed number** ("largest
  difference over the whole arc = 1.2e-14") rather than claiming it in prose.
  A checkable claim should be checked on every render.
* `density_lap_mode` splits the arc into **one trace per day**, labelled, so
  the passes can be told apart. A test asserts the split loses no points.
* The evaluation grid doubles to 721 points when the recording exceeds a
  period, so two overlapping passes are each drawn smoothly.

**The averaging that distorted the profile.** Folding a >24 h recording onto
one dial averaged the columns sharing a clock time — 08:00 on day 1 and day 2
became one number. Under extended wakefulness those two differ by the whole
homeostatic rise, so the average is a level that occurred on **neither day**;
a test asserts exactly that (the folded value lies strictly between the two
real ones and equals neither). `density_profile_mode` now offers:

* **spiral** *(new default whenever the recording spans more than a period)* —
  every column in time order as one continuous lapping path. Nothing averaged,
  nothing dropped: this is the 2-D plot in polar coordinates, and the gap
  between successive passes over a clock hour *is* the trend.
* **per_day** — one ring per day, directly comparable.
* **average** — the old fold, now saying how many clock times were averaged and
  how many observations went into the worst one.

The spiral must not be sorted or closed: sorting folds day 2 back onto day 1,
which is the averaging the mode exists to avoid, and closing joins the last
observation to the first across a gap that was never measured. It therefore
routes through `fck_open_path`/`fck_band_path` rather than
`fck_close_ring`/`fck_band_ring`, and a test asserts the clock hours are
deliberately non-monotone. The old `density_avg_days` checkbox still selects
the old behaviour, also pinned by a test.

Also: the night gradient was drawn at full strength and dominated the data.
It is now at 0.55 opacity — a background that competes with the data is a
failed background.

**4.29 Bound-hit fits are included by default, and the report says which
bound.** *(behaviour change + new table)*

4.24 / audit 2.3 excluded every fit pinned to a parameter bound from the
population summaries. That was the wrong shape of answer twice over: it
silently shrank the Circaflex sample from 1305 to 367, and it reported one
number ("938 on a bound") that could not distinguish a badly chosen constraint
from a sample full of odd subjects.

Now: **`harmonic_include_boundary`, default ON.** Non-converged fits are still
always excluded — there is no solution to average — but a bound-hit fit is a
real converged value at the edge of the feasible region, and whether to average
it is the analyst's call. The report states which choice is in force and what
it costs either way.

`fck_bounds_hit()` reports bounds **per parameter and per side**, on a relative
tolerance so a ceiling of 110 and a floor of 0.5 are judged alike;
`fck_bounds_summary()` rolls the per-subject lists into the two tables the
report prints — how often each individual bound was hit, and how many fits hit
1, 2, 3+ — plus a listing of every fit on more than one bound. Two pinned
parameters usually means a ridge: they trade off along a flat direction of the
likelihood, so neither is separately identified. `bounds_hit` and
`n_bounds_hit` also go into the parameter CSV, so the caveat survives export.

On the real data this sharpens the earlier finding considerably. It is not
"938 fits hit a bound" but:

| bound | n | % |
|---|---|---|
| `tau (upper)` | 851 | 65.2 |
| `tau (lower)` | 87 | 6.7 |

A_sat never, and **no fit on two bounds at once**. τ's ceiling there is
`t_max * 5` = 110 h — a numerical guard, not a physiological claim — so two
thirds of the sample want a τ the constraint will not give them, which is the
identifiability problem stated as a fact rather than inferred from an SD.

**4.30 Correcting the record on the sheet join.** Earlier notes and messages
said "the app joins the two sheets by row position". **That was wrong.** The
app has had an Excel sheet picker since the port (`output$excel_sheet_selector`,
ported from CIRCAREG; WaPaa always took the first sheet) and reads exactly
**one** sheet. It never opens a second and never joins.

The row-position join, and therefore the 654/410/181/59 split, came from
`tests/real_data_run.R` — my harness, not the app. The harness now defaults to
a single sheet like the app, takes `--sheet=` and `--group=`, and only performs
the cross-sheet join under an explicit `--join-sheets`, where it prints both
joins and their 38% disagreement instead of quietly picking one.

For this file there is no join to worry about: `slaperigheid` carries
`Leeftijd (in jaren)` and `Geslacht` alongside the 16 time columns, so a group
variable can be chosen from the same sheet. (Note that column holds some
non-numeric entries — `17j`, `v` — which will need cleaning before banding.)

**4.31 The time-origin toggle appeared to work backwards.** *(bug fix)*
Reported: with `Time origin = first observation` the fit plot started at 00:00
for a recording that began at 08:00, while `midnight` looked right.

The toggle was correct; the **display** was not, and only in one direction.
`first_observation` re-anchors the model, so `mod$time_vec` holds 0, 1, … 31
for a recording running 08:00 → 15:00 two days later. The axis labeller read
that as clock time and put t = 0 at midnight — eight hours early. Under
`midnight` the shift is zero and model time *is* clock time, which is why only
one of the two settings looked wrong and the whole thing read as inverted.

The fix is one addition at every **display** boundary and never inside a fit:
`fck_model_to_clock()`. Applied to the fit plot's tick text (positions stay on
the model axis, only the labels shift), its hover text, the midnight rules —
which are clock events, so on a shifted axis they no longer sit at multiples of
24 — and the polar tab, where the angle a point is drawn at is a clock time and
the lap split counts clock days. Without the shift the polar dial placed every
point eight hours early and split the laps at 08:00 instead of midnight.

The report now also prints the model axis separately from the clock axis when
they differ, so the two are never confused again, and says explicitly that the
plots convert back before labelling.

`tests/testthat/test-time-origin.R` pins both halves on the reporter's own
19-point axis: the model keeps its origin (elapsed time, and therefore the
trend, is unchanged by the shift), and everything shown to a reader converts
back — the first observation reads 08:00 under **either** setting, which is the
invariant the bug broke.

While writing those tests: `server/08_helpers_cosinor.R` used `%||%` without
defining it, relying on `00_state.R` having been sourced first. Sourced alone
by a test or a script it failed inside a helper that looked unrelated. It now
defines `%||%` when absent, so the module stands on its own.

**4.32 The model summary is collapsible.** It grew a great deal across the
audit — fit outcomes, bound tables, commonality, conditioning, per-group
detail — and all of it earns its place when checking a model and none of it
when working past one. The box is now `collapsible`, its text scrolls inside a
height set by a slider (200–2000 px) rather than pushing the rest of the tab
down, and a button saves the whole thing as a text file. The panel and the file
call one function, `.print_harmonic_summary()`, so the saved copy cannot drift
from what is on screen.

**4.33 Acrophases were reported in model coordinates as if they were clock
times.** *(bug fix — reporting layer only)*

The same class of error as 4.31, in the other half of the app. The harmonic
coefficients are fitted against the model's axis, which under
`time_origin = "first_observation"` is elapsed hours from the first
observation. An acrophase of 19.18 there means *19.18 h after 08:00* — 03:11
the next morning — and the reporting layer printed the elapsed number as a time
of day. The fitted curves never leave model coordinates and were right
throughout; only what was written next to them was wrong, and only when the
origin is non-zero.

One helper does the conversion now, `fck_acrophase_clock()` /
`fck_acrophase_label()`, and no site does its own arithmetic. The reporter's
four hand-checked values are the regression targets and all reproduce:

| model-elapsed | clock |
|---|---|
| pooled H1 19.18 h | 03:11 |
| group 1 H1 19.89 h | 03:53 |
| group 2 H1 18.93 h | 02:56 |
| H2 7.81 h | 15:49 **and** 03:49 |

**Harmonic h has h maxima per day** and all are reported. H2 repeats every 12 h,
so a single clock time is half the answer and which half matters depends on the
hypothesis.

**The H1 acrophase is not the peak of the fitted curve.** The curve also carries
the homeostatic trend and the higher harmonics, so its maximum sits elsewhere —
and it is the curve's maximum a reader sees on the plot, which is how a correct
acrophase can look wrong against its own figure. `fck_curve_peak_clock()`
reports the complete curve's maximum and minimum separately, per group as well
as pooled, and flags a maximum sitting at the edge of the observed window as the
boundary value it is rather than a peak the data contain.

**What is invariant, and therefore untouched.** A constant rotation cannot
change a dispersion or a difference. Circular SD, resultant length, the Rayleigh
Z, Watson–Williams, the Bingham/Hotelling test and every between-group contrast
are identical before and after, and tests assert it. The pairwise module is
unchanged for the same reason — a difference of two acrophases is origin-free —
but its column labels now say `elapsed h; differences are origin-free` so no one
reads a group mean off that table as a time of day.

Storage stays in model coordinates deliberately: converting at the source would
put clock times into the circular statistics, and the conversion belongs at the
display boundary where it can be applied once and audited.

**4.34 Unlabelled subjects are excluded from group analyses.** *(bug fix +
design reversal)* Reported as `Error in if: missing value where TRUE/FALSE
needed` from the group-comparison plot, with one subject missing `AGEcategory`.

4.24 / audit 1.5 carried unlabelled subjects as an `UNASSIGNED` pseudo-group so
they could not vanish silently. That was the wrong remedy: their group column is
`NA`, so `params$group == "__UNASSIGNED__"` matched nothing, the angle vector
came back empty, `circular_mean()` of an empty vector is `atan2(NaN, NaN)` =
`NaN`, and `if (NaN < 0)` is an error rather than `FALSE`.

A label-less group of one is not a group — it has no circular mean, and every
comparison built on it is undefined. They are now **excluded** from every group
analysis and counted in the audit, so the accounting still closes
(grouped + unlabelled + too-small = fitted) without poisoning the arithmetic.
Note the distinction: excluding *unlabelled* subjects is not excluding *small*
groups — a labelled group of one has a perfectly good mean direction, just an
imprecise one, and the `n >= 3` guard handles that separately. A defensive guard
also stays in the plot, since a group can still be emptied by an upstream filter.

**4.35 beta_0 is a coefficient, not a starting level.** "Intercept (β₀, at
t = 0)" invited reading the constant as the value the response started at. It is
not: at the first observation the harmonics are generally non-zero, and with a
saturating trend anchored there `S(t_min) = 0` exactly, so the two differ by the
harmonic sum. Now reported as three named quantities, pooled and per group:

* **Constant term (β₀)** — the fitted coefficient.
* **MESOR (rhythm-adjusted mean over the window)** — where the data sit once the
  homeostatic rise is counted.
* **Predicted value at the first observation** — the starting level, with its
  clock time.

All three are also stored per subject (`mesor`, `mesor_adj`, `value_at_start`)
and offered in the pairwise-comparison picker, since a group can rank
differently on each.

**4.36 Polar layout and collapsible boxes.** The acrophase polar plot's title
sat on the dial — a polar trace fills its plotting area edge to edge, so a
centred title with no reserved space lands on the 11–13 o'clock labels. The
subtitle moved to a paper-anchored annotation and the polar domain was pulled to
`y = c(0, 0.88)` to leave room, rather than shrinking the font until it stopped
colliding. The polar-density radial scale now renders horizontally
(`tickangle = 0`) and its spoke is a control rather than a constant, because
wherever it is put it will sometimes cross the data. And 49 result boxes across
13 tabs are `collapsible`.

**4.37 One group palette, keyed by name.** *(new module + fix)*
`server/02b_helpers_palette.R`. Seven palettes were in use — `rainbow()` in the
pairwise boxplots, `red/blue/green` in fPCA, Set1 in FoSR, a Brewer set in
clustering, firebrick/steelblue in the harmonic fit, a pastel set in its overlay,
and `FCK_DENSITY_COLORS` in the polar density — so the same age band was a
different colour in every figure.

The subtler half: every one indexed by **position**. Filter a group out, or sort
levels differently in one tab, and the survivors are repainted, so "the blue
group" means different things in two figures on the same screen.
`fck_group_colors()` keys on the level **name** through a stable sorted
reference, so a group keeps its colour under any filtering or reordering — the
"colour follows the entity, never its rank" rule.

The palette itself was re-stepped, not just centralised. Run through the dataviz
validator in `--pairs all` mode — the right mode here, since an overlay compares
every pair — the old one **FAILED at four groups**: amber `#eda100` against
orange `#eb6834`, normal-vision ΔE 13.7 against a floor of 15. That is exactly
the reporter's case. Slot 4 became magenta `#b5309b`, slot 5 a saturated brown
`#8a5a12`:

    light, --pairs all, 5 slots: ALL CHECKS PASS
      worst CVD ΔE 9.2 (green/orange, deutan) · worst normal-vision ΔE 19.1
      contrast WARN on the green -> relief is the legend and tables these
      figures already carry

*(Superseded by 4.41: the reporter later supplied their own palette. The
name-keyed lookup, the dash-past-the-prefix rule and everything else in this
section still hold — only the hexes and the depth of the validated prefix
changed.)*

**Five is the honest maximum.** No sixth all-pairs-separable hue exists inside
the lightness band — every candidate failed the normal-vision floor. Past five,
dash carries identity alongside hue and the figures say so. Also `rainbow()` is
gone: a rainbow ramp used as a categorical palette implies an order that is not
there and is not CVD-separable.

**4.38 Group comparisons report the MESOR, not R².** R² is a fit-quality number,
not a parameter: it says how well each subject's curve was described, not what
the group's rhythm is. Replaced by the MESOR, which sits naturally beside β₀ —
β₀ is the fitted constant, the MESOR is the rhythm-adjusted mean, and a group can
rank differently on the two.

**4.39 Warping parameters get the same ANOVA as the scores.** Registration splits
a curve into **where** it happens and **how big** it is. The component scores are
the amplitude half, computed on the registered curves; the warping functions are
the phase half, and stopping at the scores discards exactly the half that answers
"do the groups differ in timing". `fck_warp_params()` extracts the per-curve
shift and the warp amplitude (RMS distance of h(t) from the identity — defined
for any method, including nonlinear warps where no single shift exists), and they
run through the same omnibus/post-hoc machinery as **their own multiplicity
family**. Pooling them with the components would let one question borrow the
other's correction. A warp that moved nothing is reported as having nothing to
test, rather than an ANOVA on a constant.

**4.40 Smaller fixes.** The component-scores plot was capped at three PCs with a
red/blue/green palette; the cap is now a slider and the palette the shared one,
with `legendgroup` so a click hides both traces of a component. The acrophase
polar plot matches the density plot's height (540 px) so the two dials are
comparable. Three fPCA plots got a reserved top margin — a centred title with
plotly's default margin lands on the topmost trace.

Two bugs found while testing this work, both mine and both in the new code:
`fck_group_dashes()` took its reference order from the *input* order rather than
the sorted order the colours use, so a group's dash changed when levels arrived
differently even though its colour did not; and `fck_warp_params()` assigned
columns onto a zero-row `data.frame`, which throws rather than growing it.

**4.41 The palette is the reporter's "Alternating Light/Dark Pairs" set.** Four
hue families — blue `#00b0be`/`#8fd7d7`, pink `#f45f74`/`#ff8ca1`, green
`#98c127`/`#bdd373`, orange `#ffb255`/`#ffcd8e` — each with a medium and a light
member, replacing the five re-stepped hues of 4.37.

**Ordering is not the order supplied.** Taken as listed, slots 1 and 2 would be
the light and the medium *blue*: two series OKLab ΔE 14.9 apart, which read as
one group split in two rather than as two groups. A palette of this shape is
built to be used one row at a time, so the **mediums take slots 1–4 and the
lights slots 5–8**. Slot *k* and slot *k+4* are then the two members of one
family, and the common 2-, 3- and 4-group figures — including the reporter's
four age bands — get the four maximally separated hues.

**What the validator says, honestly.** In `--pairs all` light mode the four
mediums FAIL two hard checks: `#ffb255` sits at L 0.820, above the 0.77
lightness ceiling, and `#ffb255` against `#98c127` is ΔE **3.2** under
protanopia — two of the four age bands are the same colour for roughly 1 in 12
male readers. At eight slots it additionally fails the chroma floor (`#8fd7d7`
C 0.072 and `#ffcd8e` C 0.098 read as gray) and every light/medium pair within
a family falls below the 15 ΔE normal-vision floor (7.2 for the two oranges).
That last one is not a defect: those pairs are *meant* to read as related. No
subset of the eight is all-pairs separable at five slots, and no fifth hue
outside the set clears all-pairs against the four mediums — the reporter's four
families crowd three warm hues between 15° and 125°, where 4.37's palette had
magenta and a low-lightness brown for room.

**The palette ships as supplied, with the correction one constant away.**
`FCK_PALETTE_VARIANT` is `"supplied"` — the eight hexes verbatim, because that
is what was asked for and it is the reporter's call. Setting it to `"cvd_safe"`
keeps the same four families with the orange deepened (`#ffb255` → `#d38400`,
*identical* OKLCH hue angle 68.6, L 0.82 → 0.68) and the green nudged
(`#98c127` → `#9ec055`), which clears every hard check: lightness PASS, chroma
PASS, normal-vision PASS at worst 15.2, CVD at the 8.0 target boundary. Blue and
pink are untouched in both.

**Roles that are not groups moved out of the palette.** `'red'`, `'firebrick'`
and `'steelblue'` were carrying population means, single series and reference
lines across the harmonic and fPCA tabs. Red is now a problem specifically: slot
2 is a pink, so a red population-mean line reads as a fifth group.
`FCK_EMPHASIS` (dark neutral, the summary line), `FCK_NEUTRAL` (identity lines,
the confidence ellipse), `FCK_SERIES1`/`FCK_SERIES2` (the ungrouped series) and
`FCK_ALERT` (status *text*, not a series) are deliberately outside the
categorical ramp. The harmonic-component cycle
`c("green","orange","purple","brown",…)` — unvalidated CSS names — became
`fck_component_colors()`, drawn from the same ramp so one screen never mixes two
colour vocabularies.

**4.42 Three figure-layout defects.**

*The acrophase dial stacked three things at 12 o'clock.* The radial axis was
pinned at `angle = 90` — the top of the dial, once `rotation = 90` and a
clockwise direction are applied — so its "Amplitude" title and its tick numbers
were drawn straight through the topmost angular label and up into the subtitle.
The 12 angular ticks sit at 30° steps, so the radial axis moves to **45°**,
exactly midway between two of them. The subtitle became a second line of the
*title* rather than a paper-anchored annotation: at a fixed paper *y* it had no
idea where the dial ended, and plotly draws angular labels *outside* the polar
domain. The legend moved from plotly's default right edge — where a long group
name eats into the dial — to underneath it, as on the density tab. Verified by
rendering the before and after and looking at them, not by reasoning about the
geometry.

*Both polar figures reset on every keystroke.* plotly discards zoom, pan and
legend state when Shiny re-renders, so no parameter could be nudged while
watching the same view. A constant `uirevision` fixes it. The density legend
gets its **own** revision keyed on which rings exist, so hiding one group
persists while you turn the bandwidth knob but does not silently carry over to a
different set of rings when the grouping variable changes.

*The "Components to show" slider promised more than the PCA had computed.* It
ran to 10 regardless of `nharm`, so asking for 6 drew 3 and the title read "3 of
3". A display control must not outrun the analysis: the slider's ceiling is
retuned to what the PCA actually retained each time it runs, and when the two
still disagree the figure names the remedy instead of quietly drawing fewer
curves. The extraction default went 3 → 5, which is numerically inert —
functional PCA components are nested, so PC1–PC3 are bit-identical whether
`nharm` is 3 or 5 — and simply stops the common case being capped at three.

A bug found while testing this work, mine and in the test harness rather than
the app: `tests/real_data_run.R` wrote `g_pos`/`g_id` into its `.rds` output
unconditionally, but both only exist on the `--join-sheets` path.

**4.43 P0 corrections from the external statistical review.** An independent
review of the whole app raised 21 findings. I re-checked each against the code
and, where the claim was about behaviour rather than syntax, by running it. 17
confirmed, one wrong, two overstated, and one defect the review missed that was
mine. What follows is what changed and what it was measured to be worth.

*P0.2 &mdash; "Automatic smoothing (REML)" did no smoothing.* The control set
`lambda = 0` under a comment I wrote claiming `smooth.basis()` optimises
internally. It does not: it smooths with the lambda it is given and reports that
lambda's GCV score. With `nb` capped at `min(n_basis, n_time)` a few lines above,
auto mode was a square system &mdash; on 16 time points, 16 effective df, zero
residual df, data reproduced to 9e-16. Against a known truth it was 3.2x worse
than the GCV choice, and it is the default, upstream of everything.
`fck_auto_lambda()` (server/04_helpers_fd.R) now runs a real GCV search with the
same estimator that performs the smoothing: coarse log grid to find the basin,
then `optimize()` to refine, minimising mean GCV across subjects because one
lambda for the sample is what makes curves comparable downstream. On a 25-subject
fixture: df 16.00 -> 4.17, RMSE against truth 1.145 -> 0.332.

*P0.3 &mdash; the zero-amplitude test rejected true nulls three times too often.*
My own audit fix. For `exp_sat` the reduced model froze tau at its full-model
value and refit only the MESOR and A_sat by OLS. I chose that for numerical
robustness and documented the shortcut, but never worked out its cost. Under the
null tau is free, so SSE0 was not minimised, the numerator was inflated and F was
biased upwards. Measured over 3,000 simulated nulls on the Circaflex time grid:
**14.6% rejection at nominal 5%, 4.5% at nominal 1%**. `fck_reduced_exp_sat_sse()`
refits tau with `nls.lm` under the same bounds, and returns NA rather than a
number if the refit fails &mdash; a wrong F is worse than an absent one. After:
**3.6% at nominal 5%, 0.8% at 1%**. On the real data this moved "significant
rhythms" from **67.5% to 59.7%**: 102 of 1305 subjects were being called rhythmic
and are not. Any earlier report of that figure should be corrected.

*P0.4 &mdash; the repeated-measures F used a residual containing the effect.*
`Y_centered <- Y - rowMeans(Y)` removes the subject margin only, so
`sum(Y_centered^2)` was SS_error + SS_visit exactly. On a 12x4 fixture the app
returned F = 8.49 where `stats::aov` gives 37.26 &mdash; 22.8% of the correct
value, so the test was badly *conservative* and would miss real effects. Both the
observed and the permutation statistic now remove both margins and add the grand
mean back, and `df_within` is `(n-1)(k-1)` rather than `n(k-1)`. The fixed form
reproduces `stats::aov` to 1e-8.

*P0.5 &mdash; the amplitude bound was not an amplitude bound.* `b_cos` and `b_sin`
were each boxed at +/- amplitude_max, and a box of half-width A admits
`sqrt(A^2+A^2) = A*sqrt(2)`: ask for 10, get up to 14.14. The constraint region is
a disc and a box cannot express it. Refitting inside the inscribed square would
satisfy the bound but would also forbid legitimate amplitudes at off-axis phases,
so it is applied *only* when the converged fit actually violates the promise;
fits inside the disc keep the full box. `amp_bound_action` records which
happened. `amplitude_min` was threaded through six signatures and never applied;
it is now *reported* (`amplitude_below_min`) and still not imposed, because
forcing a flat subject up to a floor would invent a rhythm.

*P0.6 &mdash; degenerate fits.* A constant trajectory gave R&sup2; = 0/0 = NaN; a
perfect fit gave `log(0)`, log-likelihood `Inf` and **AIC = -Inf**, so a flat
subject won every model comparison it entered. `fck_r_squared()` returns NA when
SST = 0 (undefined, not 1 and not 0) and `fck_gaussian_loglik()` floors sigma^2
relative to the data scale. Separately, a subject with exactly one observed point
reached `approx()`, which raises; that error escaped to the module-level handler
whose fallback reverts the **whole run** to raw data for every subject. One and
two-point subjects are now handled explicitly and counted as failures so they are
reported.

*P0.8 &mdash; time warping was partly random.* `linear_shift_alignment()` added
`sin(pi*t) * runif(1, -0.03, 0.03)` to an estimated transformation under the
comment "Add slight S-curve for visualization", on a module with no `set.seed`:
measured run-to-run spread was up to 5.9x the shift being estimated. The measured
lag was multiplied by 0.1 and again by 0.5, so the warp carried 5% of what it had
measured. The landmark branch computed landmarks, discarded them, returned a warp
made entirely of random numbers, and set `registered_curves[,i] <- curves[,i]`
&mdash; it did not warp anything &mdash; while the UI announced "Time-warped PCA
completed!". The quadratic family's derivative is negative for alpha > 1 and the
default search range was c(0.5, 2), so at alpha = 2 clipping collapsed 12 of 24
grid points onto time 0.

  The shift warp is now `h(t) = t - s` with s used as measured and limited to a
  quarter of the domain; the landmark branch builds the monotone piecewise-linear
  map carrying each curve's own landmarks onto the reference ones and actually
  applies it; the quadratic family's search range is restricted to where the map
  is a bijection. `tests/warping_test.R` executes the estimator: identical output
  across runs, all warps strictly increasing, identity in gives identity out
  exactly, and a known phase shift is recovered with r = 0.998 and slope 0.901
  &mdash; against 0.045 under the old attenuation. `allow_dilation`,
  `dilation_range` and `symmetric_warp` were read from the UI, passed as
  arguments and referenced by no function body; they are labelled "not
  implemented" rather than silently accepted.

*P1.1 &mdash; permutation p-values could be exactly zero.* `mean(perm >= obs)`
reports p = 0 for a Monte Carlo test of finitely many draws. All five sites now
use `(1 + #{T* >= T}) / (1 + B)`, and the default permutation count goes from 200
(resolution 0.005, with an FDR correction on top of it) to 5,000.

**What the review got wrong.** `1 - deviance/null.deviance` labelled "McFadden
Pseudo-R&sup2;" is *correct* for a binomial GLM on 0/1 data &mdash; the saturated
log-likelihood is zero, so deviance = -2 logLik and the two are the same quantity
(verified identical to 10 decimal places). It diverges only for a proportion
response (0.652 vs 0.274), which is the separate defect the review also correctly
identifies; fixing that removes this one. The FoSR bootstrap is percentile-CI
inversion &mdash; a recognised approximate test with poor small-sample coverage,
not the meaningless number "not a valid bootstrap hypothesis test" implies; the
comment claiming it is "a proper bootstrap test" is still indefensible and the
missing multiplicity correction across time is the more serious half. The
formula-injection surface is two FoSR sites, not "several": the harmonic
`as.formula(paste(...))` calls use internal names the app controls.

**Still open** (P1/P2, and none of it changes a number already fixed above): the
mgcv -> fda lambda transfer on the diagnostics tab, the L2 statistic not
integrating over time, FoSR QR fitting and multiplicity, SoFR outcome handling
and apparent-performance labelling, clustering metrics computed in mixed
geometries, `install.packages()` at startup, and the export path writing `aov()`
where the app runs permutation + FDR &mdash; that last one means "Export analysis
code" still does not always mean "code that reproduces this result".

**4.44 P1 corrections, plus the one defect the second review found.** A second
independent review of the P0 work confirmed all eleven fixed items and raised one
finding neither review had caught. What follows is P1 as scheduled, plus that.

*P1.a &mdash; the RM permutation was drawn per time point, not per curve (NEW).*
`sample(1:n_visits)` sat INSIDE `for(t in 1:n_time)`, so a subject could be read
as condition order 2-1-3 at one time point and 3-2-1 at the next. In a functional
design the exchangeable unit is the whole trajectory.

  Being precise about what this did and did not invalidate, because the review
  did not say and it matters: at a single time point, permuting condition labels
  within a subject IS the correct null, so every pointwise p-value the app has
  reported was marginally valid. What the per-time-point draw destroyed was the
  JOINT null across time, which is what a global max-F or L2 statistic would need
  &mdash; and the RM branch sets `p_value_L2 <- NA` and computes no global test,
  so no reported number was wrong because of this. It is a real defect that
  blocked a correct global test and wasted `n_perm x n_time x n_subjects` draws;
  it is not a retraction. One relabelling is now drawn per subject per replicate
  and applied wherever that subject appears, which needs the subject indices
  behind each complete-case matrix (`Y_rows`), since that set differs by t.

*P1.b &mdash; the app claimed to run rmfanova and never did.* Four guessed
signatures in a tryCatch cascade &mdash; `rmfanova(curves, id=, visit=)`, the same
call positionally, then `exists()` probes for `rm.fanova` and `rmfANOVA` &mdash;
none of which is the package's API (it takes a list of condition-specific
matrices). The package branch could never succeed, so the app always ran its own
procedure while the module was named for the package, and refused to start
without a package it never called. Guessing at an API in a cascade is not a
fallback; it is a way of not knowing which estimator produced your numbers. The
speculative calls and the hard dependency are removed and the method is named:
pointwise repeated-measures ANOVA with within-subject, curve-wise permutation.

*P1.3 &mdash; the L2 statistic was a vector norm.* `sqrt(sum(v^2))` over the
evaluation grid scales with grid density: the same function on 50 versus 200
points gave 4.95 and 9.97. `fck_l2_norm()` uses trapezoidal weights, so the same
function gives 0.7071 on both grids &mdash; the analytic value of the L2 norm of
sin(2*pi*t) on [0,1]. Uneven grids are now weighted by their own spacing rather
than silently over-weighting dense regions. Seven sites.

*P1.2 &mdash; intervals labelled for what they are.* The fANOVA "confidence
bands" are pointwise percentile intervals, not simultaneous functional bands, and
were built from 100 bootstrap replicates (a 2.5% quantile from the 2nd-3rd order
statistic). Now 2,000 replicates, labelled "95% pointwise CI", with the UI saying
that a region where the interval excludes zero is not a family-wise claim.

*P1.6 &mdash; clustering compared quantities in different spaces.* Cluster centres
came from `values$fd_obj`, always original-scale, while the members came from
`data_matrix`, which may have been standardised. WCSS then subtracted a raw-units
centre from standardised observations. A fixture shows this drives
`between_ss = total_ss - within_ss` NEGATIVE. Objective centres now come from the
same matrix as the points; original-scale means are kept as `cluster_means_raw`
for display. Note the second review credited me with having fixed the
restart-selection WCSS &mdash; `git log` shows `60_clustering.R` untouched since
the merge, so that path was always consistent; only the final centres were wrong.

*P1.4 &mdash; FoSR.* `solve(crossprod(X))` squares the condition number and gives
no warning on a rank-deficient design; replaced with a QR that reports the rank
and names the collinear columns. The formula was built by pasting uploaded column
names into text and re-parsing; `reformulate()` takes them as data. The bootstrap
p-values (`2 * mean(boot_betas <= 0)`, commented "This is a proper bootstrap
test") resample around the fitted ALTERNATIVE, so the distribution is centred near
beta-hat: that is percentile-CI inversion, a usable interval and a poor test,
labelled the other way round. The bootstrap now yields the SE and interval only;
inference is the studentised statistic, FDR-adjusted across time &mdash; there had
been no multiplicity correction at all on a curve of ~100 pointwise tests. In the
GAM branch, coefficient curves were written only for NUMERIC predictors, so a
factor was fitted, appeared in the results, and displayed a row of zeros as
though estimated; factors now emit one contrast curve per non-reference level.

*P1.5 &mdash; SoFR.* A three-level factor response became 0, 1, 2 despite a UI
warning &mdash; a warning the code ignores is not a safeguard; now a hard error.
Bare proportions were passed to `binomial()` with no denominator, so 4/5 and
800/1000 were treated identically; now refused with the three ways out named.
`lf()` was called without `argvals`, assuming an even grid on an app that has an
explicit uneven-time option; the normalised real positions are now passed, bound
through the formula's environment rather than `pfr`'s data list (argvals has
length n_time, every data column has length n). ROC, AUC and accuracy are
labelled apparent, because they are evaluated on the fitting sample.

  On "McFadden Pseudo-R2": both reviews called this wrong. It was right on the
  branch it was written for. For a binomial GLM on genuine 0/1 data the saturated
  log-likelihood is zero, so deviance = -2*logLik and `1 - dev/null.dev` IS
  McFadden &mdash; identical to 10 decimal places on a test case. It diverges only
  for a proportion response (0.652 vs 0.274), which P1.5 now rejects outright. The
  display says "Deviance explained", which is right either way.

*P1.7 &mdash; the lambda transfer is gone.* The "Use Diagnostic Results" button
took an `mgcv::gam` smoothing parameter and handed it to `fda::fdPar()`. Those
multiply different penalties on different bases; the number transferred cleanly,
its meaning did not. The CV had a second problem: it scored lambda for predicting
a held-out subject from the population mean, which is not the problem the
smoothing module solves. Since P0.2 the app has a GCV search on the production
smoother, so the transfer is unnecessary as well as wrong. The button now runs
that search and is no longer gated on the REML/CV panels having been run; those
remain as diagnostics about smoothness and no longer set the estimator's lambda.

Also `aes_string()` -> `aes(.data[[ ]])` (2 sites), and `rmfanova` dropped from
the optional-package list since nothing calls it.

**Still open after P1** (all P2, none of it changes a number): exported code does
not call the same kernels the app runs &mdash; the fANOVA export still writes
`aov()` where the app runs permutation + FDR, which is now further from the truth
than before P0.4; `install.packages()` at startup; no `renv.lock`; Roxygen
contracts on the statistical kernels; vectorising the permutation loops now that
the default count is 5,000. The SoFR changes are the one part of this round that
is NOT runtime-verified: `refund` will not install for this container's R version,
so `lf(argvals=)` and the response guards are checked by inspection and by unit
tests on the surrounding logic, not by executing a fit.

**4.45 P2 engineering, and the rmfanova question answered with a simulation.**

*P2.1 &mdash; exported code now calls the app's own estimator.* The between-subjects
fANOVA export wrote out `summary(aov(curves_eval[, t] ~ group_labels))` and took
the parametric p-value, while the app runs a permutation test with FDR; the
repeated-measures export emitted a commented sketch of an rmfanova call matching
no version of that API. After the P0.4 residual-SS fix the exported `aov()` was
closer to the app's OLD behaviour than its current one. The harmonic export had
already solved this by deparsing the live function object into the script, so the
same pattern is applied: `deparse(perform_functional_anova)` and
`deparse(perform_rm_fanova)`. Deparsed code cannot drift from the implementation
because it IS the implementation.

*P2.2 &mdash; the app no longer installs packages while starting.* `app.R` called
`install.packages()` for every missing required package and, inside a tryCatch,
for every missing optional one. That made the estimator whatever CRAN held on the
day; it could not work on the deployed, locked-down or offline machines where it
mattered; and it hid a five-second fix inside a half-installed dependency tree.
It now checks and stops with the exact command to run. `tools/renv_bootstrap.R`
records a lockfile &mdash; deliberately NOT committed from the audit container,
because a lockfile is supposed to record the library an analysis actually ran
against and four of the app's packages are not installed here, so any lockfile
written from this machine would be a fiction.

  Exported scripts now also stamp `R.version.string` and the versions of fda,
  mgcv, refund, fda.usc, minpack.lm and the rest, and CHECK them on re-run. A
  script that names no versions is re-runnable, not reproducible.

*P2.3 &mdash; `server/01b_kernel_contracts.R`.* Both reviews said the statistical
contracts live only inside 5,800-line files. Nine kernels now have a recorded
contract: orientation, time units, estimand, estimator, assumptions,
missing-data policy, bounds, what the intervals mean, degenerate cases and
failure modes. Orientation is stated for every one because the app mixes
subjects x time and time x subjects and that is the easiest silent transpose in
the codebase. Extracting the kernels into a package is the right long-term move
and is not this change; moving eight estimators out of the reactive files at once
is how a working app stops working.

*P2.5 &mdash; the permutation loop was not runnable at the new default.* Four-deep
interpreted `for(perm) for(t) for(group) for(curve)`, which at the P1.1 default of
5,000 permutations, 100 evaluation points and the Circaflex sample is ~650 million
R-level iterations. Replaced with the matrix form: SSB is a weighted row-sum over
group means, SSW the row-sums of squared deviations from each curve's permuted
group mean. Verified identical to the loop (max difference 4e-14, floating-point
only) and ~8x faster; 5,000 permutations goes from ~47 s to ~6 s on a
100 x 400 problem.

**4.46 Is rmfanova worth implementing? Yes, with a curated default. Here is the
evidence.**

The app's repeated-measures procedure is pointwise: it says WHERE the conditions
differ and computes no overall p-value at all (`p_value_L2 <- NA`). That is a real
gap, because "do they differ at all" is usually the first question. rmfanova
(Kurylo & Smaga 2023) answers it, with three statistics (Cn, Dn, En) under five
resampling schemes plus pairwise comparisons.

*It installs, with one obstacle worth knowing.* rmfanova's DESCRIPTION lists
`refund` in Imports, but no `refund::` call appears anywhere in its source (it
uses MASS and parallel). So an install can fail for a dependency nothing needs.

*The fifteen outputs are not interchangeable.* Measured over 400 simulated nulls
and 400 alternatives (n=15, p=20, l=3, subject random intercepts, iid noise, a
Gaussian bump on one condition), at nominal 5%:

| method | type-I | power | | method | type-I | power |
|---|---|---|---|---|---|---|
| Cn_P1 | 3.8% | 100% | | Dn_B3 | 0.0% | 99.5% |
| Cn_P2 | **0.0%** | **0.0%** | | En_P1 | 6.0% | 99.8% |
| Cn_B1 | 0.2% | 100% | | En_P2 | **17.5%** | 100% |
| Cn_B2 | **0.0%** | **0.0%** | | En_B1 | 1.5% | 96.8% |
| Cn_B3 | 0.5% | 99.8% | | En_B2 | **14.2%** | 100% |
| Dn_P1 | 3.5% | 100% | | En_B3 | 4.0% | 99.2% |
| Dn_P2 | 0.0% | 76.0% | | | | |
| Dn_B1 | 0.0% | 98.2% | | | | |
| Dn_B2 | 0.0% | 73.0% | | | | |

Two reject true nulls at roughly three times the nominal rate. Two have no power
whatsoever. Handing an analyst all fifteen in a GUI is an invitation to choose
one. This is ONE data-generating scenario and not a verdict on the package &mdash;
the authors characterise these tests across many more, and a different covariance
structure could reorder them &mdash; but it is enough to justify a curated default.

*How it is wired in.* `server/09b_helpers_rmfanova.R`. Optional, off unless asked
for. It builds the list of condition matrices the documented API wants, which
needs a COMPLETE BALANCED design: every subject in every condition exactly once,
in the same row order in every matrix. The app's own pointwise procedure takes
complete cases per time point and is more forgiving, so subjects missing a
condition are dropped here and both the count and the names are reported rather
than absorbed. `FCK_RMFANOVA_DEFAULT` shows Cn_P1, Dn_P1 and En_P1; the other
twelve are in the results object and the report names which are inflated and
which have no power. `tests/rmfanova_calibration.R` re-runs the simulation and
asserts that the curated set is still the calibrated one.

**4.47 Scalar-on-function regression removed.** At the reporter's request: they
never used it, and it was the one module whose P1 fixes could not be verified by
running them, because `refund` will not install for R 4.3.3 in the audit
container. A tab that is never used, cannot be tested here, and had four
confirmed defects going into P1 is not worth carrying &mdash; every future audit
pass has to reason about it, and every future reader has to decide whether to
trust it.

Removed: `server/71_sofr.R` (702 lines), `ui/71_sofr.R` (85 lines), the sidebar
entry and tab object in `app.R`, `values$sofr_model` and its reset in
`00_state.R`, the SoFR line in the session summary, section 12 of the exported
script (which called `refund::pfr`), the `library(refund)` line the export
emitted, the `fit_sofr` kernel contract, and the SoFR fixture and section
assertion in `tests/codegen_test.R`.

Two consequences worth recording. **`refund` is no longer a dependency of the
app at all** &mdash; nothing else called it, so the package that could not be
installed here is simply gone from the requirement list. And the app's demo
generator still produces its `binary_outcome` column: it was created for the
SoFR logistic branch, but it is a perfectly good two-level grouping variable for
the fANOVA, clustering and cosinor group comparisons, so removing it would have
made the sample data worse for no reason. The comment now says why it exists.

The P1.5 tests went with the module. What replaces them is a guard that the
removal is complete &mdash; no orphan file, no menu entry pointing at a missing
tab, no `sofr_model` in state, no `refund::pfr` in the export &mdash; so it
cannot creep back half-wired. The `aes_string` guard, which had pointed at
`71_sofr.R`, now sweeps every remaining server file instead.

After removal: 12 tabs, 155 outputs, 30 server files, 44 files parsing. 1,080
testthat assertions and 10 standalone suites pass, and the Circaflex figures are
identical to P2 (59.7% significant rhythms, Rayleigh Z 549.9), confirming nothing
downstream depended on it.

**4.48 P3: eight items from a third reviewer, of which seven held.** The list
was handed over after the SoFR removal. I checked each against the code before
touching anything, on the standard set earlier in this audit &mdash; the last
review had 17 of 21 claims hold, one wrong, two overstated. This time seven of
eight held, one is refused, and acting on them turned up two further defects
neither the reviewer nor I had seen.

*P3.1 &mdash; the n-basis GCV sweep still swept at lambda = 0. Confirmed.* P0.2
removed `lam <- if(input$smooth_method == "auto") 0 else ...` from
`server/20_smoothing.R` and from the "suggest a lambda" observer. It survived
in the GCV-vs-n-basis sweep in `server/30_diagnostics.R` &mdash; the one place
in the app that calls itself a diagnostic. lambda = 0 in `fda` is the
UNPENALISED fit, not "chosen automatically", so in auto mode that tab drew the
GCV curve of a model the app would never fit, and its "optimal n_basis" was the
optimum for an unpenalised spline. The sweep now re-selects lambda by GCV at
every candidate basis size, using `fck_auto_lambda` &mdash; the same search the
app runs &mdash; on a 12-point bracket and a smaller subject cap, since it runs
once per basis size. The plot and the notification say which regime produced
the curve. The finding this exposes is worth stating: once lambda is selected,
the GCV curve is much flatter, so n_basis stops mattering above a floor.
`test-p3-corrections.R` checks that numerically (the penalised sweep's relative
range is smaller, and its GCV is never worse).

*P3.2 &mdash; FoSR should use analytical p-values and keep the bootstrap for the
interval. Confirmed.* P1.4c had already removed a bootstrap "test" that
inverted a percentile interval. What it put in its place was still mismatched:
`beta / SE_bootstrap` referred to `t_{n-p}`. A t reference is what you get when
the denominator is the analytical sigma-hat, whose scaled square is chi-square
on n-p df and independent of beta-hat; a bootstrap SE is neither, it is a Monte
Carlo estimate of the same quantity carrying its own noise of order
1/sqrt(2B) &mdash; about 5% at B = 200 &mdash; and it made a reported p-value
depend on the seed. It also bought nothing: this is a RESIDUAL bootstrap,
resampling rows of the residual matrix from the fitted homoscedastic model, so
it estimates exactly what `xtx_inv[j,j] * sigma2` estimates, with extra
variance. (A case/pairs bootstrap is the one that buys heteroscedasticity
robustness. That is a different estimator, it is not what this control has ever
run, and it is not substituted here.) The test is now always the analytical
one, FDR-adjusted across time; the bootstrap supplies the percentile interval
and its SE is kept as a misspecification diagnostic &mdash; the summary prints
the SE_bootstrap / SE_analytic range and warns when it leaves [0.8, 1.25]. The
interval and the test can now disagree at the margin, deliberately, and both
readouts say so. A test checks the p-values against `lm()` term by term and
checks they no longer move with the seed while the interval still does.

*P3.3 &mdash; the repeated-measures effect size should be partial eta-squared.
Confirmed, and it was worse than "should".* The RM branch computed the correct
repeated-measures sums of squares at every time point, used them for F, and
**threw them away**; the effect size was then rebuilt from a different,
between-subjects decomposition &mdash; `SSB / SST` over all curves, with no
subject margin anywhere in it, and over a different case set from the test
(every curve, rather than the complete cases per time point the F used). So a
repeated-measures F was reported next to a classical eta-squared, which carries
the between-subject variance the F test has already partialled out and is
therefore biased DOWN, arbitrarily far down as subject spread grows. The RM
sums of squares are now kept and partial eta-squared formed from them,
`SS_condition / (SS_condition + SS_error)`; a test verifies the identity
`df_c*F / (df_c*F + df_r)` and cross-checks F against `stats::aov`'s own
`Error(s/c)` decomposition, and shows the old number to be more than five times
too small on a design with large subject spread. Both branches now carry an
`eta_squared_type` label, the axis and title name which one is on screen, and
the classical value is still reported in the summary so earlier output can be
reconciled &mdash; explicitly marked as not to be reported.

**This changes a number in the repeated-measures tab.** The RM global L2
statistic did too: it was built from the same between-subjects SSB, and is now
the integrated condition sum of squares from the RM decomposition. It is
descriptive only (`p_value_L2` is NA for this design; the global test is
rmfanova's, behind `run_global_test`), but it is not the number earlier
versions printed.

*P3.4 &mdash; the FoSR export should be the same implementation as the GUI.
Confirmed.* Sections 7 and 10 already emitted the app's estimators with
`deparse()` of the live function object. Section 11 was a hand-written
reconstruction, and it had drifted: it emitted
`as.formula(paste('~', paste(fosr_predictors, ...)))` and
`solve(crossprod(X))` &mdash; precisely the two constructions P1.4b and P1.4a
removed from the GUI &mdash; and no standard errors, no p-values and no FDR
adjustment at all. So "Export analysis code" silently reintroduced a fixed
defect and could not reproduce any of the inference on screen. The estimator
now lives in `server/07_helpers_fosr.R` as `fck_fit_fosr_ols()`, is called by
the GUI, and is written into the script verbatim. There is no second
transcription left to drift.

*P3.5 &mdash; the exported kernels should be standalone. Confirmed, and larger
than described.* `deparse()` guarantees the exported code IS the
implementation. That guarantee stopped at the top-level function: the emitted
bodies call the app's own helpers and the script defined none of them. Measured
with `codetools::findGlobals` on the live objects, the closure is
`perform_functional_anova` &rarr; `fck_l2_norm`; `perform_rm_fanova` &rarr;
`fck_l2_norm`, `fck_rmfanova_global`, `fck_rm_design`; `fit_cosinor` &rarr;
**eleven** helpers via `fit_cosinor_nonlinear`. Worse, `perform_rm_fanova`
carried `input$rm_global_test`, `showNotification()`, `withProgress()` and
`incProgress()` in its body, so the exported script for a repeated-measures
design failed at `object 'input' not found` the moment the permutation loop
started &mdash; it could not run at all. And `covariates`, `subject_id` and
`rm_factor` were indexed by the script but never assigned anywhere in it (the
group vector was emitted truncated to ten entries with a literal `", ..."`
pasted in, which does not parse).

Fixed by: making the Shiny pieces arguments with inert defaults (`progress`,
`notify`, `run_global_test`), the app passing the real ones; and an
`emit_kernel()` that computes the dependency closure **at export time from the
live function objects**, post-order, so adding a helper call to any kernel
automatically adds its definition to the script. A hand-maintained list is the
thing that went stale in the first place. One subtlety worth recording: the
permutation loop used to sit inside `withProgress({...})`, which evaluates its
expression in the CALLER's frame; replacing that with `local({...})` would have
silently discarded every write to `F_stat_perm`. It is a bare loop now.

*P3.5b, found while doing the above &mdash; the export still called an
unpenalised fit "REML".* Two sites in `server/90_export.R` emitted
`fdPar(basis, Lfdobj = 2, lambda = 0)` under `# Smoothing method: Automatic
(REML optimization)` and `# lambda = 0 triggers automatic optimization`. This
is the same false claim P0.2 removed from the app and P3.1 removed from the
diagnostics tab, surviving in the exported script &mdash; so an auto-mode
export reproduced neither the app's smoothing nor any automatic selection. Both
sites now emit the GCV-selected lambda the app actually used, with
`fck_auto_lambda` written out so the search can be re-run on new data.

*P3.6 &mdash; commit an `renv.lock`. REFUSED, for the second time, and this
time with the refusal made executable.* A lockfile is a record of a library
that exists. `renv::snapshot()` writes the version and hash of each package as
installed on the machine it runs on, and this container has neither `renv` nor
`shinyWidgets`, `fda.usc` or `reticulate`. A lockfile written here would either
omit them &mdash; in which case `renv::restore()` silently does not install
them and the clustering and cosinor tabs fail at first use on a machine whose
owner was told the environment was pinned &mdash; or carry invented versions
and hashes, which is worse. A lockfile that does not describe a working library
is not a weaker guarantee than none; it is a false one, and this app is being
corrected precisely because it made claims of that shape. What I did instead:
`tools/renv_bootstrap.R` now reads the required/optional package lists out of
`app.R` rather than keeping a second copy that can go stale, prints the
installed version of every one, **refuses to snapshot** if any REQUIRED package
is missing, and warns which optional ones will be absent from the lock. Run on
this container it refuses, naming `shinyWidgets` &mdash; which is the point.

*P3.7 &mdash; stale SoFR and refund references in the documentation. Confirmed.*
My 4.47 removal sweep used `grep --include="*.R"` and never looked at Markdown.
`README.md` still listed scalar-on-function regression in the source-app table,
`refund` in the package table, SoFR as tab 9, and `71_sofr.R` in the file tree,
and told the reader `smoke_test.R` "runs without `fda` or `refund` installed".
All corrected, tab numbering closed up, and the test list brought up to date.
`PORTING_NOTES.md` still mentions SoFR throughout and should: it is the
chronological record of the removal. A test enforces exactly that split.

*P3.8 &mdash; a test that runs the exported script in a clean session and
compares numerically. Confirmed, and the most valuable of the eight.*
`tests/codegen_test.R` checked that the generated script PARSES. That is a low
bar: a script full of undefined symbols parses perfectly, which is how P3.5
went unnoticed. `tests/export_roundtrip_test.R` generates the script for both
designs, RUNS each in a separate `Rscript --vanilla` process on the same data
the app analysed, then re-runs every estimator from the definitions that script
supplied under the app's seeds and compares element by element. It agrees to
1e-13 on 37 quantities &mdash; the GCV lambda, the smoothed curves, F(t),
eta-squared(t), the L2 statistics, both permutation p-value vectors, twelve
cosinor fields across twelve subjects and their fitted values, and every FoSR
coefficient, SE, raw and adjusted p-value and bootstrap CI bound.

The first time it ran, it failed. The export's own cosinor reporting block read
`f$amplitude` and `f$acrophase`; `fit_cosinor()` returns those as `amplitudes`
and `acrophases`, vectors with one entry per harmonic, so both were `NULL` and
the script died on "arguments imply differing number of rows: 1, 0" at the
first subject. That defect had been in the export since the merge and no
parse-level test could ever have caught it. Fixed, along with the assumption
that every fit carries a MESOR.

After P3: 12 tabs, 31 server files, all parsing. **1,157 testthat assertions
(0 failed, 0 skipped) and 10 standalone suites pass, including the new
round-trip.** The Circaflex dataset is no longer present in this container, so
the 59.7%-significant-rhythms and Rayleigh Z figures were NOT re-measured; they
are unaffected by construction, since `08_helpers_cosinor.R` and
`72_harmonic.R` are untouched by this round. The repeated-measures partial
eta-squared and RM L2 statistic ARE new numbers and should be re-run wherever
the old ones were reported.


**4.49 P4: a third review, mostly on time warping.** Eleven findings. Nine hold
as stated, one is overstated, and one is a correction to an audit comment of
mine that was simply false. Checked each against the code before touching
anything.

*P4.1 &mdash; the parametric warp families. Confirmed, and the worst statistical
defect left in the app.* Three of the four families could not express the
identity, so a curve needing no registration was deformed anyway.

* **exponential** h(t) = (e^{at}&minus;1)/(e^a&minus;1) has its identity at the
  LIMIT a &rarr; 0. The code special-cased `abs(alpha - 1) < 0.001` and returned
  `t` there. At a = 1 the family gives (e^t&minus;1)/(e&minus;1), which differs
  from t by 0.123 at its worst &mdash; so the guard put a **discontinuity of
  that size in the middle of the default search interval [0.5, 2]**, and
  `optimize()` can converge onto it. A curve reported as "alpha = 1" had been
  left unwarped by accident while its neighbours were warped by a real map. The
  genuinely singular point, a = 0 (0/0), had no guard at all &mdash; it was
  simply outside the range the UI allowed.
* **quadratic** identity at a = 0. My own P0.8 clamp forced a &ge; 0.05, copied
  from the power family where it is correct, and the UI default gave [0.5, 1].
* **logistic** identity as steepness &rarr; 0, excluded, and 0/0 at exactly 0.
* **power** identity at a = 1 &mdash; the only family whose identity the range
  contained.

Each family now declares its identity and the open interval on which it is a
strictly increasing bijection of [0,1]; the user's range is clamped to that
interval and then *widened to contain the identity*. The UI slider starts at
&minus;5 instead of 0.1 and each family's identity is named in its label.

**Measured**, on twelve curves generated to need no registration at all (max
|h(t) &minus; t|, in units of the domain; &times;24 for hours):

| family | old range | new range | old distortion | new |
|---|---|---|---|---|
| power | [0.50, 2.00] | [0.50, 2.00] | 0.036 | 0.036 |
| exponential | [0.50, 2.00] | [0.00, 2.00] | 0.062 | 0.018 |
| quadratic | [0.50, 1.00] | [0.00, 1.00] | **0.125 (3.0 h)** | 0.018 |
| logistic | [0.50, 2.00] | [0.00, 2.00] | 0.008 | 0.007 |

For exponential and quadratic the fitted parameter came back as exactly 0.500
for every curve &mdash; the optimiser pinned at the boundary, which is the
cleanest possible evidence that the right answer lay outside the range.
`tests/warp_family_test.R` now checks endpoints, strict monotonicity, range,
exact identity and continuity in the parameter, on a 121-point grid per family.

*P4.2 &mdash; the shift warp is not endpoint-anchored. Confirmed; my comment was
wrong.* h(t) = t &minus; s maps [0,1] onto [&minus;s, 1&minus;s]. It is a
translation, which is what shift registration IS, and translations do not fix
the endpoints. The arithmetic was right and the label was wrong. What was also
wrong: a shift estimated by CIRCULAR cross-correlation was applied with
`approx(rule = 2)`, i.e. constant extrapolation &mdash; so the region the
circular estimate said should wrap round from the other end was filled with a
repeat of the endpoint value instead. On 24-hour data with a 2.4 h shift that is
a tenth of the cycle replaced by a constant. Periodic shifts now wrap;
non-periodic ones still clamp, and the extrapolated fraction is reported. The
README says plainly what each of the three registration methods is.

*P4.3 &mdash; reformulate() does not protect uploaded column names. Confirmed,
and this is a false claim I wrote in a FIX.* The P1.4b note said reformulate()
"takes the names as data and quotes them". It does not; it pastes its termlabels
and parses them, and R's own documentation says the labels must be syntactically
valid names or already backquoted. Measured: `reformulate("a b")` errors;
`reformulate("Age (years)")` silently returns `~Age(years)`, a **function
call**; `reformulate('I(cat("PWNED"))')` returns a formula holding a live call.
So P1.4b swapped one text-pasting route for another and asserted a safety
property that was not there. The GAM branch never even pretended &mdash; it
pasted names into formula text directly. Both branches now fit on internal
column names `x1..xp`, which cannot be anything but names, with the user's
labels restored afterwards as data. A test fits three hostile column names,
including `I(stop("executed"))`, and checks the coefficients match `lm()` on
renamed columns.

*P4.4 &mdash; the QR pivot. OVERSTATED as a live bug, guarded anyway.*
`chol2inv(qr.R(qrX))` is (X'X)^-1 only for an unpivoted decomposition. But
LINPACK's `dqrdc2` only cycles a column to the back when its reduced norm falls
below the rank tolerance, and it decrements the rank when it does &mdash; so
rank == ncol(X) implies pivot == 1:p, and the existing rank check has already
stopped otherwise. Measured at three conditioning levels: whenever pivoting
would have mattered, the rank test fired first. That is an undocumented
invariant of one code path, not a guarantee, and it costs one comparison to stop
relying on it, so the pivot is now checked and undone if present.

*P4.5 &mdash; n == p divides by zero degrees of freedom. Confirmed.* Such a
design is full rank, so the rank check passes, and then sigma2 = 0/0, SE = Inf
and `pt(df = 0)` = NaN, with nothing the user could connect to the cause.
Measured on a 3&times;3 design. Now refused with a message naming the counts.

*P4.6 / P4.7 &mdash; two degenerate-case guards. Confirmed.* FoSR R&sup2; had no
protection for a time point where the response is constant (NaN or &minus;Inf,
then plotted and averaged); the between-subjects fANOVA F had none for zero
within-group variation (Inf or NaN, and NaN then propagates into the permutation
comparison, where `NaN >= NaN` is NA and the p-value silently comes from fewer
draws than it reports). Both now return NA, and NA is carried through to the
p-value and excluded from the significant-region count rather than folded into
"not significant". The repeated-measures branch already did this correctly.

*P4.8 &mdash; the cyclic n-basis diagnostic. Confirmed.* The sweep always built
a B-spline basis. Under cyclic smoothing the production smoother fits a Fourier
basis, so for periodic data the diagnostic answered a question about a model the
app was not fitting. It now branches on `input$is_cyclic`, as the "suggest a
lambda" observer already did, and sweeps an ODD grid &mdash; fda rounds an even
Fourier `nbasis` up, so an even grid would have scored the same model twice
(pinned in a test).

*P4.9 &mdash; app.R pointed at an renv.lock that does not exist. Confirmed.*
Pointing at a lockfile that is not there is the same category of defect as the
rest of this audit: a claim of a guarantee that is not present. The note now
says what to run and states plainly that the project is not environment-pinned
until you run it, and the README carries a section explaining why the lockfile
must be generated on the analysis machine. My P3.6 refusal stands; the dangling
pointer was a separate and real defect.

**A pattern worth recording.** Three times now a fix has reproduced the exact
string it was removing, inside the comment explaining the removal &mdash; which
makes the grep guard that proves the removal pass for the wrong reason. It
happened again here (P4.2 and P4.9) and both were caught by the tests. The rule
is written into the code now: describe the removed text, do not repeat it.

After P4: **1,193 testthat assertions (0 failed, 0 skipped) and 11 standalone
suites pass**, including the new `warp_family_test.R` and the export round-trip.
The Circaflex data is still absent from this container; the cosinor files remain
untouched by this round, so those figures are unaffected by construction.
Anything registered with the exponential, quadratic or logistic warp families
SHOULD be re-run &mdash; those results were produced with the identity outside
the search space.


**4.50 P5: a fourth review, and a regression of my own.** Nine findings, all of
which held. One of them is a bug I introduced in the previous round's fix, and
acting on the review's suggested tests turned up a tenth defect that had been
silently wrong since the merge.

*P5.1 &mdash; FoSR prediction was broken by my P4.3 fix. Confirmed; my
regression.* P4.3 made the estimator fit on internal column names `x1..xp` so
that uploaded names could never reach a parser. It returned `terms = terms(f)`,
which therefore refers to `x1..xp` &mdash; and the three prediction sites went
on building their new-data frames with the USER's column names and calling
`model.matrix(delete.response(mod$terms), data = pred_df)`. Reproduced:
`object 'x1' not found`. Not for hostile names &mdash; **for every FoSR OLS fit,
every time**. And the call sat inside `tryCatch(..., error = function(e) NULL)`,
so the user got a blank prediction curve rather than an error. A broken panel
that looks like an empty panel is how this survived a release.

The mapping and the fitted factor levels now travel with the fit, and
`fck_fosr_design()` is the single place that turns new data into a design
matrix &mdash; used by the fit, the interactive prediction, both reference
curves and the export. The error is reported rather than swallowed.

*P5.2 &mdash; the warped-PCA AIC/BIC. Confirmed; removed.* The "residual
variance" they were built on was `mean(rmse^2)` where `rmse` is the distance
between each curve and its OWN registered version &mdash; how much the
registration moved the curve, which is what registration is FOR. As a
criterion it rewarded doing nothing. There was also no likelihood: no
probability model for the observed curves under a candidate registration exists
anywhere in the module, so `-2 log L + 2k` was applied to a number that is not a
log-likelihood. And `k_params <- 2 * n_subjects` was hard-coded regardless of
method, so the penalty did not distinguish the methods it was being used to
compare. The tab told the user to select a warping method by comparing them.
All of it is gone; the panel reports dispersion instead, and says plainly that
there is no criterion here.

*P5.3 &mdash; the "EFDA variance decomposition" was not a decomposition.
Confirmed.* Three quantities computed independently, presented as a split of a
total, under the name of a published methodology. They do not sum and nothing
established orthogonality. Renamed to pre/post registration dispersion with
`G = 1 - V_post/V_pre` defined explicitly on screen, and an explicit line
saying not to report G as "variance explained by phase".

*P5.4 &mdash; the elastic phase distance was invalid for two of the three warp
types. Confirmed, and worse than described.* Measured on a 100-point grid:

| warp | phase distance as computed |
|---|---|
| shift h(t) = t &minus; s, any s | **0.000000** &mdash; h' = 1, so the metric is identically blind to translation |
| periodic shift, wrapped, any s | **0.142254**, *including s = 0* |
| power &alpha; = 1 / 1.3 / 2 | 0.000 / 0.130 / 0.340 &mdash; behaves correctly |

The periodic case is the sharper defect: an unshifted curve was reported with
the same nonzero phase distance as a quarter-cycle shift, an artefact of the
wrap discontinuity that the `warp_deriv[warp_deriv < 0] <- 0` line converted
into a plausible-looking number. The metric is now computed per geometry:
Fisher-Rao only for endpoint-preserving warps (and only after checking the warp
IS one), displacement for a translation, circular displacement for a periodic
shift.

*P5.5 / P5.6 &mdash; landmark registration. Both confirmed.* The manual branch
had no monotonicity requirement at all; crossed landmarks produced a folded
"warp" and were registered with it. And the two branches stored **inverse maps
of each other**: the manual branch used h: registered &rarr; original, the
automatic branch built h: original &rarr; registered and then applied it as
`approx(h, curve, xout = t)`, which inverts it again. The registered curves came
out plausible either way, but `warp_functions` held two different objects
depending on which branch ran, and every statistic that reads
`warp_functions - time_points` was comparing incomparable things across
methods. One validated builder now serves both, the direction is a stated
module-wide contract stamped on all three methods' output, and rejected curves
are named in a warning instead of folded.

*P5.7 &mdash; RM-fANOVA turned undefined into "no effect". Confirmed.* P4.7
fixed this in the between-subjects branch; the RM branch still did
`F_stat[is.na(F_stat)] <- 0`, entered `0` for permutation draws with no residual
variation, and finished with `p[is.na(p)] <- 1`. Entering zero is not neutral:
zero is the SMALLEST possible F, so it never exceeds the observed statistic and
biases every p-value DOWN (pinned in a test). NA now travels through.

*P5.8 / P5.9 &mdash; the GAM branch. Both confirmed.* The export still carried a
hand-written reconstruction that pasted uploaded column names into parsed
formula text (the P4.3 defect) and specified a different model from the app's
(no basis, no k, no factor main effects). It is now a kernel,
`fck_fit_fosr_gam()`, emitted verbatim like the others, with the spline
dimension as an argument rather than three hard-coded literals. And its
`beta.se`/`beta.p` were `beta_hat * 0` labelled "placeholders" &mdash; SE = 0
asserts an estimate known exactly and p = 0 asserts overwhelming significance.
NA now, with a note saying why, and the consumer that did
`if (sum(beta.se[idx, ]) > 0)` had to be fixed too: on an all-NA row that is
`if (NA)`, which is an error, not a skip.

*P5.10 &mdash; stale REML labels. Confirmed.* The smoothing readout said
"Automatic (REML)" and "Lambda: 0 (automatic REML optimization)" over a GCV
smoother, and the diagnostics tab told the user lambda = 0 uses REML. This is a
documentation defect with a real consequence in a statistical application: the
user writes the wrong method in a paper. Corrected everywhere, with the mgcv
REML panels relabelled as advisory and on a different scale.

*P5.11 &mdash; the parametric search was using an optimiser whose assumption the
objective violates. FOUND BY THE TEST THE REVIEW ASKED FOR, not by the review.*
The reviewer suggested a registration-effectiveness test: simulate phase-only
data and check that registration reduces between-curve dispersion. Written, and
it failed immediately &mdash; the power family gave **G = &minus;7.4%** on
phase-only data (registration made it worse) and deformed already-aligned curves
by **0.58 of the domain**.

The cause is not the families, which P4.1 fixed. It is `optimize()`. Golden
section plus parabolic interpolation assumes a UNIMODAL objective. The
registration SSE on a sharply peaked curve &mdash; which is what a circadian
profile is &mdash; is not: measured on an aligned sample, SSE is 0.008 at
&alpha; = 1 and about 20 everywhere else in [0.05, 6], a deep narrow well in a
wide plateau, and `optimize()` returned &alpha; = 6.000. **This got worse at
P4.1, not better**: widening the ranges to make each identity reachable was
right, but a wider interval gives `optimize()` more plateau to get lost on, so
registrations run with the narrow pre-P4.1 ranges were partly protected by luck.
Replaced with a coarse grid scan plus refinement in the winning bracket &mdash;
the pattern `fck_auto_lambda()` already used.

| | before | after |
|---|---|---|
| G, phase-only data (power) | &minus;7.4% | **+99.8%** |
| deformation of aligned curves (power) | 0.5816 | **0.0010** |

*P5.12 &mdash; amplitude leakage, and a knob I declined to add.* With the search
fixed, one case remained: on curves differing ONLY in amplitude, the logistic
family reports G = 27.8% from warps averaging 0.014 of the domain, peak heights
unchanged. That is amplitude being absorbed as phase &mdash; near a peak a 1%
move in time changes the value a lot &mdash; and it is intrinsic to
least-squares registration, not a coding error. The standard remedy is a
deviation-from-identity penalty. I implemented and calibrated one, and **it does
not work here**: across &lambda; from 0 to 0.2 it moves G by 1.4 points, because
the offending warps are already near-identity so there is nothing for the
penalty to bite on. Shipping it would have been a tuning knob that looks like a
fix. Instead the panel now WARNS on the signature &mdash; a large G from a
near-identity warp &mdash; and the README states the measured example.

After P5: 32 server files, **1,309 testthat assertions (0 failed, 0 skipped) and
12 standalone suites pass**, including the new
`tests/registration_effectiveness_test.R`. The Circaflex data is still absent
from this container; the cosinor files are untouched by this round.

**Re-run anything registered with a parametric warp.** P5.11 changes the fitted
parameter on any data whose alignment objective is multimodal, which includes
every peaked profile. Anything read off the warping panel's AIC/BIC or its
"variance explained by warping" should be discarded outright &mdash; those
numbers did not mean what the panel said.


**4.51 P6: seven bugs found by running the app.** Not a review this time --
someone opened the app on real data and it misbehaved in seven places. Two of
them are mine, from the previous two rounds, and *none* of them was visible to
any test in this repository. That is the important finding, and it is addressed
at the end.

*P6.1 &mdash; FoSR GAM: "object 'j' not found". My regression, from P5.8.* When
the GAM branch became a kernel I renamed its loop indices from `k` to `j`,
because `k` shadowed the new spline-dimension argument. One loop got the
rename in its BODY but not its HEADER: `for (k in seq_along(preds))` with
`preds[j]` inside. The GAM branch died on its first line of real work.

*P6.2 &mdash; fPCA-ANOVA: "$ operator is invalid for atomic vectors".* Two R
behaviours compounding.
`res$warping_method <- NULL` does not store NULL in a list, it **deletes the
element**; and `$` on a list falls back to **partial matching**. With
`warping_method` deleted, the only remaining name beginning with "warp" was
`warped`, so `res$warp` returned `res$warped` &mdash; a logical scalar.
`!is.null(res$warp)` was therefore TRUE, the phase-comparison section ran, and
`res$warp$k` died. It fired on every run where warping had NOT been used, i.e.
the common case. (With both `warped` and `warping_method` present the prefix is
ambiguous and `$` returns NULL, which is why the warped path was fine and the
bug looked intermittent.) Fixed twice over: the slots are always created, with
NA rather than NULL, so no dangling prefix exists; and every read of that family
uses `[["..."]]`, which is exact.

*P6.3 &mdash; every fPCA display capped at three components.* Ask for five, get
three, in the summary, the loadings plot, the effect-of-scores plot and the
slider's ceiling. Four places asked "how many components are there?" as
`length(pca_res$harmonics)`. **`harmonics` is an fd object, and an fd object is
a list of three elements**, so that expression returns 3 for every PCA ever run.
Every consumer took `min(..., 3)`. The summary had a hard-coded
`1:min(3, length(varprop))` on top of that. Now `fck_n_harmonics()`, which reads
the coefficient columns, and the summary prints every retained component with a
running total.

*P6.4 &mdash; the "Components to show" slider drew ~19 ticks labelled
1,1,1,2,2,2,2,3,3,3 on a 1-to-3 range.* `updateSliderInput()` was sent only the
new `max`; the client then recomputed tick positions from the old range. Send
min, max and step together.

*P6.5 &mdash; the FoSR observed-data plot showed one group out of four.*
`subset_idx <- 1:min(nrow(Y), 200)` &mdash; the first 200 rows in FILE ORDER.
Data exported per group arrives sorted by group, so on 654 YOUTH / 410 ADULT /
181 MIDDLE_AGE / 59 ELDERLY the plot drew 200 YOUTH curves and one legend
entry, and looked exactly like a dataset with one category in it. The cap is
there to keep the plot drawable; it must not decide WHICH curves you see. Now a
stratified draw with a floor of five per level, a fixed seed so the picture is
stable, and a subtitle saying what was sampled. It also used
`scale_color_brewer("Set1")` &mdash; a fourth palette in an app that has one;
now the app's own.

*P6.6 &mdash; six pairwise p-values, all exactly 0.009901.* Not a bug: with
B permutations the smallest attainable Monte Carlo p is 1/(B+1), and at B = 100
that is 0.009901. Six clearly-significant comparisons all land on the floor and
the table shows the identical number six times, which is indistinguishable from
a bug. The summary now states the floor, counts how many comparisons sit on it,
and says the identical values are the resolution limit rather than a tie. The
pairwise permutation default went from 200 (min 100) to 5,000 (min 500) to match
the omnibus tab, which P1.1 already raised.

*P6.7 &mdash; three defaults the wrong way round.* Data source defaulted to
smoothed, when cosinor is a regression on the observations and handles gaps
natively &mdash; smoothing first inflates R&sup2;, makes LOOCV optimistic and
biases the zero-amplitude F test anticonservatively, all of which the app's own
help text already said. Time origin defaulted to midnight, which with a
saturating trend gives an intercept that is the value at neither origin. "Show
raw data points" defaulted on, burying the fitted curve on a few hundred
subjects. All three flipped, and the server-side `%||%` fallbacks flipped with
them &mdash; a session that never touches a control must not run a different
model from the one the UI displays.

**The gap this exposed, and what was done about it.** Every test in this
repository was either static (does it parse, does the source contain X) or
numerical (does this kernel return the right number for this input). Nothing
pressed a button. So a renamed loop variable and an R partial-match both shipped
&mdash; neither is a statistical error, and both made the app unusable in a tab.

`tests/reactive_smoke_test.R` is the missing layer. It builds a dataset shaped
like the real one (sorted by group, unbalanced, four levels), drives the server
with `shiny::testServer` through fPCA, the component ANOVA and both FoSR
methods, and FORCES EVERY RELEVANT OUTPUT TO RENDER; an output that errors fails
the test. It was verified by reintroducing P6.1 and P6.3 one at a time and
confirming it fails on each, then restoring.

One thing it surfaced immediately and worth recording: several plot outputs wrap
their bodies in `tryCatch(..., error = function(e) cat(...))`, so a broken plot
renders as an empty panel and prints to a console the user never sees. That is
the same failure mode that let P5.1 survive a release. It is not fixed here, but
it is now visible: the reactive test drives those outputs, so a future error in
them shows up as text in the test log even when the app would have swallowed it.

After P6: **1,372 testthat assertions (0 failed, 0 skipped) and 13 standalone
suites pass.**


**4.52 P7: a control that had never done anything.** Reported directly: typing
9904 into the pairwise tab's "Number of permutations" produced a run at 1000.

*P7.1 &mdash; the pairwise permutation box was dead code.* The observer read

```r
n_perm_to_use <- if (!is.null(values$fanova_results$n_permutations))
                   values$fanova_results$n_permutations
                 else input$pairwise_permutations   # "Fallback"
```

under the comment "IMPORTANT: Use the same n_permutations that was used in the
omnibus FANOVA". **The tab refuses to run without an omnibus fANOVA**, so the
first branch is always taken and the second is unreachable. The control has
never had any effect; the readout reported the omnibus's count, and nothing said
the typed value had been discarded. P6.6 made it worse by adding help text
telling the user to raise it.

There is no statistical reason the two counts must match &mdash; they are
different tests, and the pairwise family needs MORE resolution, because a
correction across *m* comparisons pushes the smallest attainable adjusted
p-value to *m*/(B+1). The control is honoured now, both counts are reported so a
difference is visible, and a blank or invalid entry falls back with a warning
rather than silently.

While there: the summary now also reports the smallest attainable **adjusted**
p, and the app **refuses the run** when *m*/(B+1) exceeds alpha &mdash; at
B = 100 with six comparisons and Bonferroni that is 0.059, so nothing could be
significant however strong the effect, which previously came out as six null
results with no explanation.

*P7.2 &mdash; a sweep for the same shape found two more.* Comparing every
`*Input("id")` and `actionButton("id")` in `ui/` against every `input$id` in
`server/` turned up `max_landmarks` ("Maximum number of landmarks", enforcing
none) and `display_options` ("Show individual curves" / "Show mean curves",
read by no plot). Both are wired.

**Two checks, because neither subsumes the other.** The static sweep is now a
test &mdash; it would have caught P7.2 and it asserts the list stays empty. It
would NOT have caught P7.1, because `input$pairwise_permutations` *was*
referenced, in a branch that never ran. So `tests/reactive_smoke_test.R` also
sets the control to a distinctive value (777, against an omnibus at 200) and
asserts the result carries it. Verified by reintroducing the override and
confirming the test fails with
"the pairwise tests used B = 200, not the 777 that was asked for".

After P7: **1,387 testthat assertions (0 failed, 0 skipped) and 13 standalone
suites pass.**


**4.53 P8: a fifth review. Four findings, all of which held, and the first is
mine again.**

*P8.1 &mdash; GAM PREDICTION still used the user's column names. Confirmed; my
omission.* P4.3 made both FoSR branches fit on internal names `x1..xp`, and
P5.1 fixed the consequent prediction failure &mdash; **for the OLS branch
only**. `build_gam_pred_df()` went on keying its frame by the user's names,
against a formula that refers to `x1..xp`. Reproduced:

```
ERROR: object 'x2' not found
Warning: not all required variables have been supplied in newdata!
```

on every GAM prediction curve and both min/max reference curves. The mapping is
now a REQUIRED argument of that helper, so a caller who has not thought about
it cannot call the function at all; all three call sites pass
`mod$gam_model_names`.

The reviewer also identified exactly why it survived: *"Your new GAM test
currently tests fit, not fit + prediction."* True. `tests/reactive_smoke_test.R`
now drives the real helper, predicts, and checks the frame is keyed by the
model's own names &mdash; and that a call without the mapping is refused.

*P8.2 &mdash; the periodic-shift warp amplitude was a representation artefact.
Confirmed, and it reached the statistics.* A wrapped warp
h(t) = (t &minus; s) mod 1 is a CORRECT registration of a cycle &mdash; t = 1
and t = 0 are the same instant &mdash; but drawn on [0,1] it contains one jump,
and an RMS distance from the identity reads that jump as deformation. Measured
on a 100-point grid:

| true shift | RMS &#124;h &minus; t&#124; reported |
|---|---|
| 0.00 | **0.1000** |
| 0.05 | 0.2179 |
| 0.10 | 0.3000 |
| 0.25 | 0.4330 |

An unshifted curve carried the warping amplitude of a real 0.1 shift &mdash;
2.4 hours on a 24-hour day &mdash; and the values are inflated and non-monotone
in s. The velocity variance is worse: var(h') = 98 at every shift, zero
included. **And this was not only a display number**: `fck_warp_params()` fed
`warp_amplitude` into the group comparison on warping parameters, so a "do the
groups differ in phase" test on periodic-shift data was testing the artefact.

A translation's intensity is its displacement &mdash; the shortest CIRCULAR one
when the design is periodic. Only an endpoint-preserving warp is measured by its
distance from the identity. Velocity variance is NA for a periodic translation.
One consequence worth recording: for a shift, `warp_amplitude` is now
&#124;shift&#124;, a deterministic function of the `shift` column beside it, so
the redundant column is dropped rather than spending a second entry of the
multiplicity family on the same question.

*P8.3 &mdash; a ratio between two incomparable lambdas. Confirmed.* The
diagnostics printed `cv_lambda / reml_lambda` and concluded "REML and CV agree
well on optimal smoothing". They are penalty weights on DIFFERENT penalties over
DIFFERENT bases; the ratio is not a measure of agreement and a value near 1 does
not mean the methods concur. Worse, the report converted the mgcv REML optimum
to `-log10(lambda)` and offered it as a "Smoothing Factor" &mdash; a number the
user could type into a control where it means something else. Both removed.
(The CV lambda IS on the fda scale &mdash; the CV loop calls
`fdPar(basis, 2, lambda)` and `smooth.basis` &mdash; so its smoothing factor
remains transferable, and stays.)

*P8.4 &mdash; the CV panel's estimand was mislabelled. Confirmed.* That CV
smooths the TRAINING-GROUP MEAN and scores it against a held-out subject, so it
answers "how much should the population curve be smoothed to predict a new
person". The production smoother answers "how much should this person's own
trajectory be smoothed". Different questions; recommending the first as
"lambda for prediction tasks" invited the reader to use it for the second. The
panel is now titled *Between-subject population-curve prediction CV* and says
so explicitly.

*P8.5 &mdash; the UI was contradicting its own output. Confirmed.*
`ui/41_results.R` still announced "Warping Fit Statistics ... based on EFDA
methodology", "Variance Decomposition" and "Model Selection Criteria" directly
above panels that (since P5.2/P5.3) say there is no AIC or BIC and that the
dispersion numbers are not a decomposition. Retitled to Registration
Diagnostics / Pre-post registration dispersion / Registration summary, with the
per-subject table relabelled as what it is: R&sup2;, RMSE and MAE there compare
each curve with its OWN registered version, so they measure the magnitude of the
transformation, and doing nothing scores perfectly. The README's summary section
still said "automatic REML or manual lambda" and "warping fit statistics and
variance decomposition" while its own detailed section said the opposite; both
synchronised.

On renv.lock, the reviewer now accepts the code is honest and asks only that a
lockfile be generated before release. That is the position taken at P3.6 and it
has not changed: it must be generated on the analysis machine, and
`tools/renv_bootstrap.R` refuses to write one from a library that is missing
required packages.

After P8: **1,440 testthat assertions (0 failed, 0 skipped) and 13 standalone
suites pass.** Anything read off the warping panel for a PERIODIC linear-shift
registration should be re-run &mdash; the warp amplitudes, the velocity
variances, and any group comparison on warping parameters were computed on the
artefact.


**4.54 P9: a sixth review. Four findings, all of which held.** No new
statistical blocker; these are the last of the "the code says one thing and the
label says another" class, plus one more instance of the recurring pattern.

*P9.1 &mdash; the per-subject "Warping Amplitude Scores" table still carried the
periodic-shift artefact. Confirmed. THE THIRD COPY.* P8.2 corrected "how much
was this curve warped" in the registration statistics and in
`fck_warp_params()`. There was a third body, in the table, and it was missed.
Every linear-shift result carries `warp_functions`, so its first branch always
won and the table went on reporting a PERIODIC ZERO SHIFT as 0.1 of the domain.
Its third branch, `abs(alpha_values - 1)`, was wrong for a different reason:
since P4.1 only the power family has its identity at alpha = 1.

The reviewer's framing is the right one &mdash; *"You've repeatedly found that
duplicated statistical definitions drift apart; this is another example."* It is
the fourth time in this audit. There is now one definition,
`fck_warp_amplitude()` in `server/04_helpers_fd.R`, called by the table, the
registration statistics and the group comparison. Verified numerically: a
periodic zero shift returns exactly 0, a 0.9 shift returns 0.1 (0.9 forward is
0.1 backward), and an endpoint-preserving warp still returns its RMS distance
from the identity.

*P9.2 &mdash; the diagnostics UI still said the two lambdas should agree.
Confirmed.* P8.3 removed the ratio and the "REML and CV agree well" verdict from
the server, and left `ui/30_diagnostics.R` telling the user *"Compare REML vs
CV: Should agree on optimal range"* and *"If REML and CV agree: use their
recommended lambda"*. The app could state both positions at once. Replaced with
what each panel is for, and a plain instruction not to compare them numerically.

The reviewer also caught something I should not have needed telling: **internal
audit tags had leaked into user-facing text**. "(P5.10)" appeared three times in
the diagnostics UI, and "(P3.2)", "(P3.3)" and "(P5.9)" in output the user
reads &mdash; the last of them emitted into the generated analysis script. They
belong in source comments. A test now scans every string literal in `ui/` and
`server/` for the pattern and requires none.

*P9.3 &mdash; the CV lambda was not on the production scale after all.
Confirmed, and it makes P8.4's "you can type this in" wrong.* The CV built
`create.bspline.basis(rangeval = c(1, n_time), nbasis = nb)` on the integer
column index. That is production's basis only in the default case: under cyclic
smoothing production fits a FOURIER basis, and on real clock times it fits over
elapsed hours on an uneven grid where the roughness penalty is per hour. In
either case the CV's lambda weighted a different penalty over a different basis
&mdash; the same category of error as the mgcv/fda ratio, one level less
obvious, and I had just certified the transfer as sound.

The reviewer offered a simple option (stop claiming transferability) and a
preferred one (make the CV build the production basis). The preferred one is
taken, because it also removes the duplicated definition rather than papering
over it: `fck_smoothing_axis()` and `fck_smoothing_basis()` now serve both, and
a test checks the extracted builder reproduces all three production branches
exactly &mdash; including that a cyclic basis on real elapsed hours gets a
24-hour period rather than the length of the recording. The CV result records
which basis it used and the report prints it.

*P9.4 &mdash; an unavailable standard error was reported as 0. Confirmed.*
`safe_gam_predict()` returned `se = rep(0, ...)` on its fallback path, and the
FoSR prediction did the same when the SE vector came back the wrong length. SE =
0 is the strongest possible claim &mdash; the prediction known exactly &mdash;
returned precisely when the standard error could not be computed. Downstream it
drew a confidence band of zero width. NA now, the band is omitted, and the plot
says why: a missing band and a zero-width band look identical and mean opposite
things. Same rule as P5.9 for the GAM coefficient curves.

After P9: **1,536 testthat assertions (0 failed, 0 skipped) and 13 standalone
suites pass.** Anything read off the per-subject Warping Amplitude Scores table
for a periodic registration should be re-run.


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
* **Not run end-to-end here** *(as of the original port; largely superseded by
  4.43&ndash;4.48).* `tests/smoke_test.R` and `tests/clock_helpers_test.R` verify
  that every file parses, the whole UI renders, every tab is reachable, no
  output is defined twice and all server files register in one environment.
  They cannot verify the statistics. At the time of the port `fda`, `rmfanova`,
  `fda.usc` and `shinyWidgets` were not installable here, so nothing had been
  executed since the merge. Since then `fda`, `mgcv`, `minpack.lm`, `rmfanova`
  and `lme4` have been installed and the estimators ARE exercised &mdash; by the
  testthat suite, by `tests/export_roundtrip_test.R`, and by the calibration and
  real-data runs recorded in 4.43&ndash;4.48. `fda.usc`, `shinyWidgets`,
  `reticulate` and `circular` are still missing, so the clustering tab and the
  UI widgets that need them remain unexecuted. Run each tab once against a known
  dataset before trusting the numbers.
