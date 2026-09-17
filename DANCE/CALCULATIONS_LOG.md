# Calculations log — cosinor module overhaul

This is the log requested with the overhaul: which of the numbers the cosinor
module prints come straight from an R package, which are computed in this
repository, and which of those in-house computations are verified by a test
in `tests/`. It is a developer document, not app text.

The rule the module follows: no procedure is invented here. Everything is
either a package call, a textbook formula implemented directly, or plain
arithmetic on package output. The table below is the audit of that rule.

## 1. What the two approaches compute

### Mixed-effects approach (one model)

| Quantity | Where | Source of the calculation | Verified by |
|---|---|---|---|
| The model fit | `server/08e_helpers_trajfit.R` `dance_traj_fit()` | `lme4::lmer` / `lmerTest` (REML); optional `glmmTMB` for AR(1)/OU residuals | `tests/traj_validation_grid.R` (type-I rates), `tests/traj_framework_test.R` |
| Convergence / singular / rank-deficient status | `dance_traj_status()` | read off the `merMod` object (`lme4::isSingular`, optimizer messages, rank of X) | `tests/testthat/test-p16-corrections.R` |
| Random-effects ladder (fall-back to simpler structures) | `dance_traj_re_ladder()` | in-house *strategy* (which structures to try, in what order); every rung is still an `lmer` fit | `tests/traj_framework_test.R` |
| Cell curves and pointwise bands | `server/08g` `dance_traj_predict()` | fixed-effect design row × β, SE from `vcov(model)`; z multiplier | `tests/testthat/test-traj-components-polar.R` |
| Simultaneous (Scheffé) bands | `dance_traj_band_mult()` | Scheffé multiplier √(q·F<sub>q,ddf</sub>) — textbook, implemented here; ddf = smallest Satterthwaite df of any coefficient (`lmerTest`) | not separately simulated — **in-house choice of ddf**, stated in the band note |
| Cell (cos, sin) coefficients and their covariance | `server/08f` `dance_traj_cell_coefs()` | `emmeans::emtrends` on the c<sub>h</sub>/s<sub>h</sub> columns | `tests/traj_corrections_test.R` |
| Cell amplitude, acrophase | `dance_traj_cell_coefs()`, `dance_traj_cell_equations()` | √(a²+b²), atan2(b, a) — plain arithmetic on the emmeans output | `tests/testthat/test-traj-participant-cell.R` |
| Amplitude / acrophase intervals, pairwise amplitude & phase comparisons | `dance_traj_contrasts(method = "joint")` | **in-house**: joint multivariate-normal draws of the fixed-effect vector (β̂, V̂) and the functional evaluated per draw; percentile intervals; circular difference on the shortest arc. The delta-method alternative (`method = "delta"`) is also implemented. | `tests/testthat/test-traj-pairwise-blocks.R` (structure and orientation, not coverage). **Not simulated for coverage.** |
| Hotelling-type joint region on (a, b) per cell, and the "acrophase not identified" rule | `dance_traj_cell_coefs()` / polar dial | 2-df ellipse from the emmeans covariance; identifiability = ellipse excludes the origin (Bingham et al. 1982 logic) | `tests/testthat/test-traj-components-polar.R` |
| Omnibus tests of the design factors (whole trajectory, shape, circadian, trend, level blocks) | `server/08f` `dance_traj_marginal_test()` → `dance_traj_L_test()` | **in-house L-matrix construction** (difference contrasts on the factors in the effect, equal weights over the others — the emmeans convention), tested with `pbkrtest::KRmodcomp(m, L)` (Kenward–Roger) or `lmerTest::contest(m, L, joint = TRUE)` (Satterthwaite); glmmTMB: Wald χ² from `vcov` | `tests/traj_corrections_test.R`, `tests/traj_validation_grid.R` (type-I rate of the circadian block on the marginal path: 0.045 / 0.050 at N = 1000) |
| Estimability check of each L | `dance_traj_L_estimable()` | row space test against the model matrix — standard linear-model algebra | `tests/traj_corrections_test.R` |
| Nested fixed-effect comparisons (ML refits) | `dance_traj_refit()` | `lme4`/`glmmTMB` refit with `REML = FALSE`; likelihood-ratio via `anova()` | `tests/traj_corrections_test.R` |
| Calibration statement ("outside the simulated grid") | `dance_traj_calibration()` | rule-based comparison of the fit's configuration with the configurations the validation grid covered | `tests/testthat/test-p16-corrections.R` |
| Participant estimates (individual table, BLUP cloud on the polar plot) | `server/08g` `dance_traj_participant_table()`, `dance_traj_participant_curves()` | own cell's fixed effects + `lme4::ranef` conditional modes (subject, and curve when present) — the standard BLUP; per-curve predictions via `predict(model, re.form = NULL)` | `tests/testthat/test-traj-participant-cell.R` (planted 12 vs 5 amplitude recovered per cell; the earlier grand-mean bug is the case this test pins) |
| Marginal and conditional R² | `dance_traj_r2()` | **in-house implementation** of Nakagawa & Schielzeth (2013) with Johnson's (2014) random-slope extension (mean over observations of zᵢ'Σzᵢ), from `getME(m, "Z")`, `VarCorr`, `sigma`. Neither `performance` nor `MuMIn` is a declared dependency. | `tests/testthat/test-traj-r2.R`: equals the hand computation from the model's own Σ and design rows (random intercept, and random slopes where the diagonal-only shortcut would be wrong), and equals `performance::r2_nakagawa()` on the random-intercept model when that package is present |
| Residual SD, logLik, AIC, BIC | summary panel | `stats::sigma`, `stats::logLik`, `stats::AIC`, `stats::BIC` on the `merMod` | — |
| Curve maximum / minimum ("complete fitted curve peaks at …") | `dance_traj_curve_peaks()` | grid search (1441 points over one period) on `dance_traj_predict()` output — arithmetic, resolution 1 minute at a 24 h period | `tests/testthat/test-traj-acrophase-display.R` |
| Residuals for the diagnostics tab | `harmonic_residuals()` | `residuals(model)` (conditional residuals) | `tests/reactive_smoke_test.R` renders the panel |
| Polar density (mixed) | `server/74_polar_density.R` | kernel density of the BLUP acrophases (same density code as before, different input); the fit rings are `dance_traj_predict()` on the dial's grid, i.e. the tab-1/tab-6 curves | `tests/polar_agreement_test.R`, `tests/testthat/test-polar-laps.R` |

### Two-stage approach (one OLS cosinor per participant, then compare)

| Quantity | Where | Source of the calculation | Verified by |
|---|---|---|---|
| Per-participant fit | `server/72_harmonic.R` `fit_cosinor()` | `stats::lm` on the cosinor design (linear, log trends); `minpack.lm::nlsLM` / `nls(port)` for the saturating trend | `tests/circular_inference_test.R`, `tests/reactive_smoke_test.R` |
| Zero-amplitude test | `fit_cosinor()` | nested-model F (full vs trend-only) — textbook | `tests/testthat/test-p6-runtime-bugs.R` |
| Bingham joint region per participant | `server/08b` `dance_bingham_ci()` | Bingham et al. (1982) ellipse on (a, b) from the OLS covariance | `tests/pop_cosinor_test.R` |
| Group vector means (amplitude-weighted) | `server/08` | mean of the (a, b) vectors — Bingham's group mean | `tests/circular_inference_test.R` |
| MESOR (rhythm-adjusted mean over the window) | `dance_rhythm_adjusted_mean()` | numerical integration of β₀ + S(t) over the observed window — **in-house definition** (the constant term is not the MESOR when a trend is present) | `tests/testthat/test-p14-corrections.R` |
| Nested-model comparison (ΔAICc, Akaike weights) | `dance_model_selection()` | AICc per participant from `stats::logLik`, averaged; weights = exp(−Δ/2) normalised — textbook | `tests/reactive_smoke_test.R` |
| Free-vs-fixed-tau check | `conditioning` block | ΔAIC between two `nls` fits per participant, averaged | `tests/reactive_smoke_test.R` |
| PRIMARY group test | `server/08h` `dance_ts_manova()` | `stats::manova`, Wilks' λ (`summary(..., test = "Wilks")`); one-way ANOVA when only one coefficient | `tests/testthat/test-twostage-helpers.R` (equal to a direct `stats::manova` call) |
| Scalar component tests (constant, MESOR, trend) | `dance_ts_components()` | `dance_group_linear_test()` = one-way ANOVA (`stats::oneway.test`-equivalent, ω² effect size) | `tests/circular_inference_test.R` |
| Per-harmonic rhythm tests | `dance_ts_components()` | `dance_pop_cosinor()` = Bingham et al. (1982) population-mean cosinor (joint (a, b) MANOVA; amplitude and acrophase F tests; the amplitude-interpretability caution, checked) | `tests/pop_cosinor_test.R` (simulated null rates) |
| Acrophase test, unweighted | `dance_ts_components()` | `dance_watson_williams_test()` (Watson–Williams F) with the von Mises concentration check `dance_ww_assumption()` | `tests/circular_inference_test.R` |
| Pairwise comparisons | `dance_ts_pairwise()` | Welch's t (`stats::t.test`) with Cohen's d; Watson–Williams for acrophase; `stats::p.adjust` (Holm / Bonferroni) | `tests/testthat/test-twostage-helpers.R` (equal to a direct `t.test` call; planted equal/unequal pairs) |
| Group curves and bands (fitted-curves tab, comparison tab AND the polar-density rings) | `dance_ts_group_curves()` | line = curve of the group's mean coefficients (the pooled `pop_mean_fit` coefficients when there is no grouping); band = ± z · SE(t), SE(t) = SD across participants of their own fitted curves / √n — **in-house, and only pointwise**. It replaces the old band on the fitted-curves tab and on the polar-density rings, which scaled the whole curve by the SE of the H1 amplitude and was not an interval for anything. | `tests/testthat/test-twostage-helpers.R` (line and SE(t) checked against a direct computation), `tests/traj_tab6_test.R`, `tests/polar_agreement_test.R` |
| Group difference curve (comparison tab) | `ts_diff()` in `server/72_harmonic.R` | difference of two group lines, SEs added in quadrature — arithmetic | `tests/traj_tab6_test.R` |
| Bootstrap CIs (optional) | run observer | participant-level percentile bootstrap; circular percentile interval for acrophase | `tests/testthat/test-p21-*.R` |
| Commonality analysis (variance shares) | `pop_mean_fit$indiv_means` | Chevan & Sutherland (1991) partition from three nested R² — textbook arithmetic | `tests/testthat/test-p14-corrections.R` |

## 2. Calculations that are in-house and NOT backed by a package

These are the items a reviewer should know are this repository's own
implementation rather than a package's:

1. **Nakagawa/Johnson R²** (`dance_traj_r2`). Published formula, own code,
   checked against a hand computation and (random-intercept case) against
   `performance::r2_nakagawa()`. If `performance` becomes a dependency the
   function should be replaced by that call.
2. **Joint-draw intervals for amplitude, acrophase and their pairwise
   differences** (`dance_traj_contrasts(method = "joint")`). Standard
   parametric-bootstrap-of-the-coefficients logic; coverage has not been
   simulated. The delta-method alternative is available and agrees on the
   point estimates.
3. **The marginal L-matrix construction** for the block tests
   (`dance_traj_marginal_L`). It reproduces emmeans' equal-weight marginal
   convention by hand so that the SAME L can be handed to `KRmodcomp` and
   `contest`. The validation grid measured its type-I rate for the circadian
   block only (0.045 / 0.050 at N = 1000); the other blocks, multi-harmonic
   and trend configurations are flagged "outside the simulated grid" on
   screen, in the summary and in the report.
4. **The two-stage SE(t) band** (`dance_ts_group_curves`). A pointwise
   standard-error band from the spread of participants' fitted curves. It is
   descriptive; no simultaneous coverage is claimed and the note says so.
5. **The Scheffé band's denominator df** (smallest Satterthwaite df among the
   coefficients). Conservative by construction; not simulated.
6. **The random-effects ladder** (which simpler structures to try when the
   maximal one does not converge). A strategy, not a statistic; every fit on
   it is an ordinary `lmer` fit, and the rung reached is printed everywhere.
7. **Kenward–Roger cost estimate** on the comparison tab (the "may take on
   the order of N s" line). A heuristic from measured timings; it decides
   nothing.
8. **The curve-peak search** (grid of 1441 points per period). Arithmetic on
   the fitted curve; its resolution (1 min at 24 h) is stated.
9. **MESOR as the window mean of β₀ + S(t)** under a trend. A definition
   this app makes explicit; the constant term is printed beside it.

Nothing in the list above is a new inferential procedure. Items 1–5 are
either published formulas coded here or descriptive constructions whose
limits are stated at the point of use.

## 3. What was removed as obsolete in this overhaul

- `server/73_cosinor_pairwise.R`, `ui/73_cosinor_pairwise.R` — the
  "Cosinor: pairwise tests" module; its tests now live in
  `server/08h_helpers_twostage.R` behind the comparison tab under the
  two-stage approach.
- The legacy two-stage block of the comparison tab (four group-statistics
  variants) and the "Parameter distribution" tab.
- The data-source and time-origin controls (always raw observations, always
  t = 0 at the first observation).
- `dance_mixed_cosinor()` in `server/06_helpers_mixed.R`, its readout branch
  and its export branch: unreachable from the UI since the trajectory model
  replaced it.
- The `two_process` trend branches in `server/72_harmonic.R` (nine sites):
  no UI option produced that trend type and the helpers those branches
  called (`compute_mean_S_from_fits`, `predict_two_process_mean_curve`,
  `compute_two_process_trend_line`) did not exist anywhere in the repository
  — the branches would have errored had they ever been reached.
- The old fitted-curves band on the two-stage path and the same construction
  on the polar-density rings (whole curve scaled by the H1 amplitude SE).
- `dance_cosinor_param_label()`, `dance_cosinor_origin_label()` and the
  `hp_*` state slots read by the publication report.
- Roughly 600 duplicated lines in the fitted-curves renderer (the same
  group/pooled/component logic written four times) — now one block that
  draws what `dance_ts_group_curves()` returns.

## 4. Tests added or changed for the overhaul

- `tests/testthat/test-traj-participant-cell.R` — per-cell BLUP baseline
  (the grand-mean bug).
- `tests/testthat/test-twostage-helpers.R` (60 expectations) — the
  `dance_ts_*` kernels on planted data (frame, MANOVA, components, pairwise,
  curves and SE band).
- `tests/testthat/test-traj-cell-equations.R` (26 expectations) — the printed
  cell equation equals the drawn cell curve, in both coefficient and
  amplitude/acrophase form, and the (cos, sin) pairs equal the emmeans
  extraction.
- `tests/testthat/test-traj-r2.R` (16 expectations) — the in-house R² against
  a hand computation and against `performance`.
- `tests/testthat/test-p12-corrections.R`, `test-p13-corrections.R`,
  `test-p18-corrections.R`, `test-p6-runtime-bugs.R`,
  `test-harmonic-ui-structure.R` — retargeted from the removed module and
  controls to the approach-aware report and UI.
- `tests/testthat/test-p15-corrections.R`, `test-palette.R` — source-string
  assertions retargeted to the rewritten summary and the longer polar
  renderer.
- `tests/reactive_smoke_test.R`, `tests/mixed_design_test.R` — drive the
  two-stage approach explicitly; the mixed cosinor kernel block is gone.

Known, pre-existing and unrelated: five expectations in
`tests/testthat/test-palette.R` ("warping parameters are extracted…", "warp
amplitude works for a nonlinear warp…") fail when that file is run ALONE and
pass when the whole `testthat` directory runs in one session. They fail the
same way on the commit before this overhaul (checked in a clean worktree);
`dance_warp_params()` and its file were not touched here.

See the commit message for the sweep results at the time of the change.
