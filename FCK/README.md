# F\*CK — **F**unctional data analysis, **C**ircadian regression, **K**-means clustering

One Shiny app combining two previously separate tools:

| source app | what it did | where it is now |
|---|---|---|
| `WaPaa1_3.R` — *Functional Data Analysis Suite* | fPCA / time-warped PCA, functional ANOVA (between-subjects and repeated measures), post-hoc tests, functional clustering | the **F** and **K** tabs |
| `CIRCAREG.R` — *Functional Regression Suite* | function-on-scalar regression, harmonic (cosinor) regression, pairwise tests on circadian parameters | the **C** tabs |

Both apps started with the same four steps — load a file, pick the variables,
smooth the curves, check the smoothing — implemented twice. **In this app those
four steps exist once.** Import and smooth your curves, then run any of the
analyses on them without re-importing, re-selecting or re-smoothing, and
without any risk that the fPCA and the cosinor fit were run on differently
smoothed data.

Everything downstream of smoothing is the original code, carried across by line
range, so each analysis computes and prints exactly what it did in its own app.
See `PORTING_NOTES.md` for the full record of what was merged, what was
renamed, and the handful of deliberate behaviour changes.

## Running it

```r
shiny::runApp("FCK")          # from the repository root
```

Required: `shiny`, `shinydashboard`, `shinyWidgets`, `fda`, `mgcv`, `plotly`,
`DT`, `dplyr`, `tidyr`, `ggplot2`, `cluster`, `readxl`.

Optional — each backs exactly one feature, and the app starts and reports which
are absent rather than refusing to launch (both source apps aborted at startup
if any of these were missing):

| package | feature it backs |
|---|---|
| `rmfanova` | repeated-measures functional ANOVA |
| `fda.usc` | functional k-means clustering |
| `reticulate` | DCF (density-core-finding) clustering via Python |
| `minpack.lm` | robust exponential-saturation cosinor fits |
| `gridExtra`, `viridis` | plot export, continuous colour scales |

## The tabs

**Shared pipeline — do these once**

1. **Data Import** — CSV/TXT/TSV/Excel (with sheet selection), wide or long,
   separator auto-detected or chosen. One variable-selection step defines the
   curve columns *and* the scalar variables.
2. **Data Preprocessing/Smoothing** — B-spline or (for 24-hour data) Fourier
   basis, automatic REML or manual lambda, per-subject smoothing that
   interpolates missing values, optional clamping to a range, per-subject
   R²/RMSE/EDF/GCV, with a control for how many curves the plot draws.
   Optionally spaces the time points by their **real clock
   times** instead of the column index — tick this if your measurements are
   unevenly spaced (hourly by day, 2-hourly at night), or a long gap gets
   smoothed as though it were a short one.
   The **Missing data & filled points** panel below it answers "which of these
   values did I actually measure?": a map of every cell as observed /
   interpolated / extrapolated, a per-subject table sorted worst-first, and a
   one-subject inspector showing the curve against its real points. By default
   the curve is **not** extrapolated past a subject's first and last
   measurement — it is held flat there, because a spline carried past its data
   follows its end polynomial.
3. **Smoothing Diagnostics** — GAM REML fit, REML profile over lambda,
   k-fold cross-validation (optionally stratified by group), GCV vs n-basis
   sweep, and a "use these results" button that sets the smoothing factor.

**F — functional data analysis**

4. **fPCA/time-warped PCA Settings** and **5. Functional PCA Results** —
   components, loadings, variance explained, scores; time-warping with linear
   shift / parametric / landmark alignment, warping fit statistics and variance
   decomposition.
6. **Functional ANOVA** — between-subjects and repeated-measures designs,
   pointwise and global permutation tests, effect sizes.
7. **fANOVA: post-hoc tests** — pairwise curve comparisons with multiple-testing
   correction, difference/p-value plots, heatmap, significance timeline. A
   control at the top says **what** is being compared: by default the same
   variable, levels, curves and design the omnibus ANOVA used, or a variable and
   design chosen here — including a within-subjects (paired) comparison, which
   does not require the omnibus to have been a repeated-measures one.

**C — circadian and functional regression**

8. **Function-on-Scalar (FoSR)** — pointwise OLS (analytical SE and p-values,
   FDR-adjusted across time; optional residual bootstrap for a percentile
   interval), or smoothed OLS via `mgcv::gam`; coefficient, p-value, R² and
   residual curves.
9. **Harmonic Regression** — cosinor with 1–3 harmonics, optional linear/log/
    exponential-saturation trend, MESOR / amplitude / acrophase per subject,
    polar plots, circular statistics, group comparisons. The **Polar Density**
    tab draws a filled shape on a clock face — noon at the top, midnight at the
    bottom, night shaded — carrying the fitted cosinor curves from the
    Fitted Curves tab (same coefficients, same band), a von Mises kernel density
    of the acrophases (so a distribution straddling midnight reads as one peak,
    not two), or the signal itself averaged over the clock.
10. **Cosinor: pairwise tests** — pairwise group comparisons of any cosinor
    parameter, with corrections, effect sizes and confidence intervals.

**K — clustering**

12. **Functional Clustering** — k-means (`fda.usc::kmeans.fd`), hierarchical and
    DCF clustering, elbow and silhouette diagnostics, cluster mean curves,
    cluster-by-group composition tests, dendrograms.

13. **Data Export** — every table above as CSV, plot bundles as PDF, the
    smoothed curves in wide and long form, a reproducible R script covering
    every analysis family, and **save/restore of the whole session** as a
    single `.rds` (data, smoothing, every fitted model, and the package
    versions they were computed under).

## Data contract

Everything downstream reads these, and only these:

| object | what it holds |
|---|---|
| `values$data` | numeric matrix, subjects × time points, as selected |
| `values$smooth_data` | the same matrix after the shared smoothing step |
| `values$fd_obj` | `fda` object for the smoothed curves, on a 0–1 range |
| `values$time_labels` | original column names, in file order |
| `values$time_numeric` | column indices `1:n_time` — WaPaa's plotting x axis, *not* clock times |
| `values$time_clock` | real clock hours parsed from the column names, or `NULL` when they cannot be trusted |
| `values$covariates` | scalar variables, original types (predictors/response) |
| `values$group_variables`, `$selected_group_vars`, `$group_labels` | the same scalar variables as factors, for grouping |

## Layout

```
FCK/
  app.R                    packages, UI assembly, server assembly
  ui/                      one file per tab
    00_theme.R             shared CSS
    10_import.R  20_preprocess.R          <- hand-merged shared steps
    30_diagnostics.R … 90_export.R        <- ported tabs
  server/
    00_state.R             the shared state bus
    01_helpers_time.R      WaPaa's time helpers (used by every plot)
    02_helpers_gam.R       GAM prediction helpers
    03_helpers_clock.R     real clock times, when they can be trusted
    04_helpers_fd.R        the one rule for building an fd object
    05_helpers_missing.R   which values are measured, filled, or invented
    06_helpers_posthoc.R   what the post-hoc tests compare
    07_helpers_circular.R  circular density and the clock-face orientation
    10_import.R  20_smoothing.R           <- hand-merged shared steps
    11_import_views.R  21_smoothing_views.R  30_diagnostics.R
    40_fpca.R  50_fanova.R  60_clustering.R
    70_fosr.R  72_harmonic.R  73_cosinor_pairwise.R
    90_export.R            exports + the reproducible-code generator
    91_session.R           save / restore a whole session
  tests/                   see "Checking it still assembles" below
```

The server files are sourced into **one** environment, exactly as when each app
was a single script — so a helper defined in one file is visible to all of
them, and no analysis code had to be rewritten to be namespaced.

## Checking it still assembles

```r
Rscript tests/smoke_test.R           # from the FCK directory
Rscript tests/clock_helpers_test.R
Rscript tests/codegen_test.R
Rscript tests/export_roundtrip_test.R
Rscript tests/missing_status_test.R
Rscript tests/posthoc_source_test.R
Rscript tests/circular_density_test.R
Rscript tests/warping_test.R
Rscript tests/warp_family_test.R
Rscript tests/registration_effectiveness_test.R
Rscript tests/audit_test.R
Rscript tests/polar_agreement_test.R
Rscript -e 'testthat::test_dir("tests/testthat")'
```

`smoke_test.R` parses every file, builds and renders the whole UI, checks that
every sidebar entry reaches a uniquely-named tab, checks that no output id is
assigned twice, and registers all 31 server files under a mock session. It
needs only the UI packages, so it runs without `fda` or `fda.usc` installed.

`codegen_test.R` checks that the exported script *parses* for every combination
of results the generator might be asked about. `export_roundtrip_test.R` goes
further: it RUNS the exported script in a clean `Rscript --vanilla` session on
the same data the app analysed, then re-runs each estimator from the
definitions that script supplied and checks every coefficient, standard error,
p-value and permutation statistic against the app's to 1e-8. Parsing is a low
bar — a script full of undefined symbols parses — so the round-trip is what
actually backs the reproducibility claim. It needs `fda` and `minpack.lm`.

`clock_helpers_test.R` pins the clock-time parsing that the cosinor tab's
shared-times option and the real-time smoothing both depend on — including the
midnight unwrap and the cases that must be *refused* rather than mistaken for
hours. It needs no packages at all.

`missing_status_test.R` pins the observed / interpolated / extrapolated
classification the missing-data panel rests on. It needs no packages.

`posthoc_source_test.R` pins what the post-hoc tests compare — that they follow
the omnibus ANOVA's variable rather than the first scalar variable selected, and
that they refuse rather than guess when nothing lines up. No packages needed.

`circular_density_test.R` pins the clock-face orientation and the wrapping of
the acrophase density. No packages needed.

`codegen_test.R` drives the code export with a stub result for every analysis
family and checks that what comes out is valid R — a script that does not parse
looks like a reproducibility guarantee and is not one.

## What the registration methods are, exactly

The time-warping controls on the fPCA tab are three different things, and the
difference matters for what you can claim about the result.

**Parametric alignment** fits one of four one-parameter families —
`power` t^α, `exponential` (e^{αt}−1)/(e^α−1), `quadratic` αt²+(1−α)t, or
`logistic` (normalised sigmoid). Each is a strictly increasing map of [0,1] onto
itself with h(0)=0 and h(1)=1: a genuine reparameterisation of time.
`tests/warp_family_test.R` checks that on a fine parameter grid for every
family. **The identity is reachable in all four** — at α=1 for `power` and
α=0 for the other three — so "this curve needs no registration" is an answer
the search can give. It could not, before: three of the four families had their
identity outside the allowed range, and on curves needing no registration the
quadratic family deformed the time axis by up to 12.5% of the domain (3 hours
of a 24-hour day) with the fitted parameter pinned at the range boundary.

**Linear shift** estimates a lag by cross-correlation and applies h(t)=t−s.
That is a *translation* — the standard shift-registration model — and a
translation is **not** an endpoint-preserving diffeomorphism: it maps [0,1] onto
[−s, 1−s]. What happens at the boundary depends on the design. Tick *periodic*
and the warp wraps, which is correct for a full cycle of circadian data; leave
it unticked and the ends are filled by constant extrapolation, and the fraction
of each curve that was extrapolated is reported. Shifts are capped at a quarter
of the domain, beyond which a lag is not identified.

**Landmark alignment** builds a monotone piecewise-linear map through the
matched landmarks with the endpoints pinned. Landmarks that cross or duplicate
cannot define a time warp; those curves are **left unregistered** and named in a
warning, rather than being registered with a fold in them.

All three return the warp in one direction — `h` maps **registered time to
original time**, and registration is always
`registered <- approx(time_points, original, xout = h)` — so `h(t) - t` means
the same thing whichever method produced it.

**What the warping panel reports, and what it does not.** It reports `V_pre`
and `V_post`, the between-curve dispersion before and after registration, and
`G = 1 - V_post/V_pre`. G is the relative reduction in dispersion. It is **not**
an amplitude/phase variance decomposition — nothing here establishes that the
two are additive or orthogonal — and it is not labelled as one. There is **no
AIC or BIC**: no likelihood for the observed curves under a candidate
registration is written down anywhere in this module, so there is no criterion
to select a method with. Choose a registration method from what you know about
the data.

Two honest caveats the panel will tell you about. A Fisher–Rao phase distance is
computed only for endpoint-preserving warps; for a shift it is undefined (a
translation has `h' = 1` everywhere) and the panel reports the displacement
instead. And a large G obtained from a near-identity warp is flagged as possible
**amplitude leakage**: near a peak a 1% move in time changes the value a lot, so
least-squares registration can absorb amplitude differences as phase. Measured on
curves differing only in amplitude, the logistic family reports G = 28% from
warps averaging 0.014 of the domain, with peak heights unchanged. A
deviation-from-identity penalty does not fix this — the offending warps are
already near-identity — so the panel names the signature instead of hiding it.

Not implemented: SRVF / Fisher–Rao elastic registration, and any separation of
amplitude from phase variance beyond what these three give you. If your question
is specifically about phase, prefer the acrophase from the cosinor tab, which is
estimated rather than assumed.

## Reproducibility: pin the environment yourself

This repository does **not** ship an `renv.lock`, and the app is not
environment-pinned until you make one. That is deliberate. A lockfile records a
library that exists — `renv::snapshot()` writes the version and hash of each
package as installed on the machine it runs on — and it can only honestly be
written where the app's packages are actually installed. A lockfile listing
packages that were never there would make `renv::restore()` reproduce the gap
instead of filling it, on a machine whose owner had been told the environment
was pinned.

So, once, on the machine you analyse on:

```r
Rscript tools/renv_bootstrap.R   # from the FCK directory; writes renv.lock
```

It reads the required/optional package lists out of `app.R` rather than keeping
a second copy, prints the installed version of every one, **refuses to write a
lockfile if any required package is missing**, and warns you which optional
features will be absent from the lock. Commit the `renv.lock` it produces
alongside your analysis. Later, or elsewhere:

```r
renv::restore()
```

The exported analysis script is pinned separately and always: it records the
version of every package it used and prints a notice on re-run if any of them
has changed.

## A note on smoothing

There is exactly one place a roughness penalty is applied: **Apply Smoothing**
on the preprocessing tab. Nothing downstream re-smooths — the fPCA, fANOVA,
post-hoc and clustering tabs all reuse `values$fd_obj`, which is built from the
already-smoothed curves with λ = 0 (a change of representation, not a second
smooth).

If you never smooth, the curves are represented by an **interpolating** basis:
they pass through your data points and nothing is smoothed. The app says so
when it does this. (In WaPaa, both "Raw data (no smoothing)" and going straight
to an analysis quietly projected onto ≤ 20 basis functions instead.)
