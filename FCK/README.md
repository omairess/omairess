# F\*CK — **F**unctional data analysis, **C**ircadian regression, **K**-means clustering

One Shiny app combining two previously separate tools:

| source app | what it did | where it is now |
|---|---|---|
| `WaPaa1_3.R` — *Functional Data Analysis Suite* | fPCA / time-warped PCA, functional ANOVA (between-subjects and repeated measures), post-hoc tests, functional clustering | the **F** and **K** tabs |
| `CIRCAREG.R` — *Functional Regression Suite* | function-on-scalar regression, scalar-on-function regression, harmonic (cosinor) regression, pairwise tests on circadian parameters | the **C** tabs |

Both apps started with the same four steps — load a file, pick the variables,
smooth the curves, check the smoothing — implemented twice. **In this app those
four steps exist once.** Import and smooth your curves, then run any of the
seven analyses on them without re-importing, re-selecting or re-smoothing, and
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
| `refund` | scalar-on-function regression (`refund::pfr`) |
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
   R²/RMSE/EDF/GCV. Optionally spaces the time points by their **real clock
   times** instead of the column index — tick this if your measurements are
   unevenly spaced (hourly by day, 2-hourly at night), or a long gap gets
   smoothed as though it were a short one.
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
   correction, difference/p-value plots, heatmap, significance timeline.

**C — circadian and functional regression**

8. **Function-on-Scalar (FoSR)** — pointwise OLS with residual bootstrap, or
   smoothed OLS via `mgcv::gam`; coefficient, p-value, R² and residual curves.
9. **Scalar-on-Function (SoFR)** — `refund::pfr`, Gaussian or binomial, with
   ROC, calibration and classification metrics for binary outcomes.
10. **Harmonic Regression** — cosinor with 1–3 harmonics, optional linear/log/
    exponential-saturation trend, MESOR / amplitude / acrophase per subject,
    polar plots, circular statistics, group comparisons.
11. **Cosinor: pairwise tests** — pairwise group comparisons of any cosinor
    parameter, with corrections, effect sizes and confidence intervals.

**K — clustering**

12. **Functional Clustering** — k-means (`fda.usc::kmeans.fd`), hierarchical and
    DCF clustering, elbow and silhouette diagnostics, cluster mean curves,
    cluster-by-group composition tests, dendrograms.

13. **Data Export** — every table above as CSV, plot bundles as PDF, the
    smoothed curves in wide and long form, and a reproducible R script.

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
    10_import.R  20_smoothing.R           <- hand-merged shared steps
    11_import_views.R  21_smoothing_views.R  30_diagnostics.R
    40_fpca.R  50_fanova.R  60_clustering.R
    70_fosr.R  71_sofr.R  72_harmonic.R  73_cosinor_pairwise.R
    90_export.R
  tests/smoke_test.R       structural check (see below)
```

The server files are sourced into **one** environment, exactly as when each app
was a single script — so a helper defined in one file is visible to all of
them, and no analysis code had to be rewritten to be namespaced.

## Checking it still assembles

```r
Rscript tests/smoke_test.R           # from the FCK directory
Rscript tests/clock_helpers_test.R
```

`smoke_test.R` parses every file, builds and renders the whole UI, checks that
every sidebar entry reaches a uniquely-named tab, checks that no output id is
assigned twice, and registers all 17 server files under a mock session. It
needs only the UI packages, so it runs without `fda` or `refund` installed.

`clock_helpers_test.R` pins the clock-time parsing that the cosinor tab's
shared-times option and the real-time smoothing both depend on — including the
midnight unwrap and the cases that must be *refused* rather than mistaken for
hours. It needs no packages at all.
