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
| Omnibus tests of the design factors (whole trajectory, shape, circadian, trend blocks; intercept and MESOR) | `server/08f` `dance_traj_marginal_test()` → `dance_traj_L_test()` | **in-house L-matrix construction** (difference contrasts on the factors in the effect, equal weights over the others — the emmeans convention), tested with `pbkrtest::KRmodcomp(m, L)` (Kenward–Roger) or `lmerTest::contest(m, L, joint = TRUE)` (Satterthwaite); glmmTMB: Wald χ² from `vcov` | `tests/traj_corrections_test.R`, `tests/traj_validation_grid.R` (type-I rate of the circadian block on the marginal path: 0.045 / 0.050 at N = 1000) |
| Estimability check of each L | `dance_traj_L_estimable()` | row space test against the model matrix — standard linear-model algebra | `tests/traj_corrections_test.R` |
| Nested fixed-effect comparisons (ML refits) | `dance_traj_refit()` | `lme4`/`glmmTMB` refit with `REML = FALSE`; likelihood-ratio via `anova()` | `tests/traj_corrections_test.R` |
| Calibration statement ("outside the simulated grid") | `dance_traj_calibration()` | rule-based comparison of the fit's configuration with the configurations the validation grid covered | `tests/testthat/test-p16-corrections.R` |
| Participant estimates (individual table, BLUP cloud on the polar plot) | `server/08g` `dance_traj_participant_table()`, `dance_traj_participant_curves()` | own cell's fixed effects + `lme4::ranef` conditional modes (subject, and curve when present) — the standard BLUP; per-curve predictions via `predict(model, re.form = NULL)` | `tests/testthat/test-traj-participant-cell.R` (planted 12 vs 5 amplitude recovered per cell; the earlier grand-mean bug is the case this test pins) |
| Marginal and conditional R² | `dance_traj_r2()` | **in-house implementation** of Nakagawa & Schielzeth (2013) with Johnson's (2014) random-slope extension (mean over observations of zᵢ'Σzᵢ), from `getME(m, "Z")`, `VarCorr`, `sigma`. Neither `performance` nor `MuMIn` is a declared dependency. | `tests/testthat/test-traj-r2.R`: equals the hand computation from the model's own Σ and design rows (random intercept, and random slopes where the diagonal-only shortcut would be wrong), and equals `performance::r2_nakagawa()` on the random-intercept model when that package is present |
| Intercept, MESOR and non-periodic change per cell, and their pairwise differences | `server/08f` `dance_traj_cell_functionals()`, `dance_traj_level_summary()`; `server/08d` `dance_traj_trend_window()` | **in-house closed forms** for the trend basis over the window [t0, t_max], T = t_max − t0: linear mean T/2, end T; log mean ((T+1)log(T+1) − T)/T, end log(T+1); saturating mean 1 − (τ/T)(1 − e^(−T/τ)), end 1 − e^(−T/τ). MESOR = β₀ + Σ c_j·mean_j, change = Σ c_j·end_j — exact linear functions of β, so the per-cell Wald interval is from `vcov` and the tests are the same L-matrix contrasts as every block (`KRmodcomp` / `contest`). With one trend coefficient per cell the change test IS the coefficient's test | `tests/testthat/test-traj-levels.R` (closed forms against numerical integration; MESOR against the trapezoid integral of the fitted midline; change test identical to the trend test), `tests/traj_levels_app_test.R` |
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
| Non-periodic change over the recording, per participant | `fit` loop in `server/72_harmonic.R` (`trend_change`), `server/08h` `dance_ts_trend_change()` | S(t_last) − S(t_first) through `dance_trend_value()`, from the participant's own trend parameters — for the saturating trend A and τ together | `tests/testthat/test-twostage-helpers.R`, `tests/traj_levels_app_test.R` |
| Group intercept, MESOR and change (summary, report, tab 6) | `group_fits`, `pop_mean_fit`, `server/08h` `dance_ts_level_summary()` | arithmetic means of the per-participant values (the population-mean cosinor's convention), t intervals on n − 1 df; pairwise Welch's t as for every scalar | `tests/testthat/test-twostage-helpers.R`, `tests/traj_levels_app_test.R` |
| Nested-model comparison (ΔAICc, Akaike weights) | `dance_model_selection()` | AICc per participant from `stats::logLik`, averaged; weights = exp(−Δ/2) normalised — textbook | `tests/reactive_smoke_test.R` |
| Free-vs-fixed-tau check | `conditioning` block | ΔAIC between two `nls` fits per participant, averaged | `tests/reactive_smoke_test.R` |
| PRIMARY group test | `server/08h` `dance_ts_manova()` | `stats::manova`, Wilks' λ (`summary(..., test = "Wilks")`); one-way ANOVA when only one coefficient | `tests/testthat/test-twostage-helpers.R` (equal to a direct `stats::manova` call) |
| Scalar component tests (intercept, MESOR, change, trend parameters) | `dance_ts_components()` | `dance_group_linear_test()` = one-way ANOVA (`stats::oneway.test`-equivalent, ω² effect size) | `tests/circular_inference_test.R` |
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
   this app makes explicit; the intercept is printed beside it. Mixed: the
   closed forms of section 5; two-stage: numerical integration per
   participant, then the mean.
10. **The non-periodic change over the recording** (S at the last observation
   minus S at the first). Plain arithmetic on the fitted trend; in the mixed
   model a fixed multiple of the trend coefficient (section 5).

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

## 5. Intercept, MESOR and the non-periodic change

Three level quantities per cell (mixed) or per group (two-stage), shown under
*Estimated trajectories* on the comparison tab and offered in the pairwise
comparisons:

- **Intercept** — the midline (intercept + trend, without the harmonics) at
  the first observation, where every trend column is zero. β₀.
- **MESOR** — that midline averaged over the recording. Without a trend it is
  the intercept, and only the intercept is offered.
- **Change** — the trend's value at the last observation: the non-periodic
  change over the recording, the right-hand end of the "change from origin"
  view.

There is one MESOR per curve, not one per harmonic: the harmonics are
deviations about one shared midline, and with a trend each cycle of a
harmonic sits on a different stretch of that midline, so a per-harmonic MESOR
would need a definition this module does not invent.

**Mixed-effects.** With the trend basis anchored at t0 and τ fixed, each trend
column has a closed-form window mean and end value (section 1), so the three
quantities are rows of the map from β and are tested by the existing L-matrix
code. The model estimates ONE trend coefficient per cell (β for linear and
log, the rise A for the saturating trend), so the change is that coefficient
times a constant: its test is the coefficient's test, expressed in the
response's units. τ is held at the value set under Model Specification.

**Two-stage.** Each participant's fit estimates both A and τ, so the
participant's change, A·[(1 − e^(−(t_last−t_off)/τ)) − (1 − e^(−(t_first−t_off)/τ))],
summarises both in one number; A and τ remain separately comparable.

Corrected in this change:

- The mixed model's τ was read from the two-stage "free-vs-fixed" field,
  which is hidden under the mixed approach, so the fit silently held τ at
  that field's default (18 h). It now has its own visible input
  (`harmonic_traj_tau`), restored with a saved session.
- The two-stage pooled and group MESOR integrated the trend over the
  clock-hour window (e.g. 8 to 33) instead of the model window (0 to 25)
  whenever the first observation was not at midnight. On the simulated
  four-group design in `tests/traj_levels_app_test.R` the ADULT group read
  43.72 instead of 40.92.
- The group MESOR is now the mean of the participants' own MESORs, the
  quantity the group tests compare. The MESOR of the group's mean parameters
  is kept as `rhythm_adjusted_mean_of_mean_coefs`; for the saturating trend
  it differs (38.79 vs 40.92 on the same data) because the MESOR is
  nonlinear in τ.
- The Export tab's cosinor summary was a separate text that printed β₀ as
  "MESOR" and ignored the mixed approach; it now writes the on-screen
  summary.


## 6. tau: profiled for all cells, or estimated per cell with A and A/tau

**The trajectory model's tau** (mixed approach, saturating trend) is either held
at the value set under Model Specification or estimated once for all cells by
profile likelihood (`dance_traj_profile_tau`, `server/08e`):

- ML refits of the conditionally linear model on a 24-point log grid of tau
  (span/50 to 3 x span), **all on one random structure** -- the one that
  converges at the middle of the grid. Each grid point used to walk the
  random-effects ladder on its own, so a point fitted on a simpler structure
  had its log-likelihood compared with its neighbours' from a different model.
- The maximum is refined between the best point's neighbours by a
  one-dimensional search on log tau (`stats::optimize`, ~1% resolution); the
  refined point is added to the profile.
- The 95% profile interval is where the log-likelihood falls 1.92 below its
  maximum (chi-square(1)), interpolated linearly in log tau between grid points;
  a side that never falls that far is reported as open.
- A flat profile (log-likelihood range < 1.92) means tau is not identified: the
  helper returns no fit at its argmax, and the app holds tau at the value on
  screen and says so.
- The grid fits use `lme4::lmer` directly (`light = TRUE`): same likelihood as
  `lmerTest::lmer`, half the time; the final REML fit at the estimate is an
  ordinary `lmerTest` fit, and its tests are conditional on the estimated tau.

**Per cell: A, tau and A/tau** (`server/08i_helpers_trajsat.R`), on request.
A nonlinear mixed model, `nlme::nlme` (Pinheiro & Bates 2000), maximum
likelihood:

    y = b0[c] + A[c] (1 - exp(-t / exp(lt[c]))) + sum_h a_h[c] cos + b_h[c] sin
        + participant random effects (diagonal; ladder: level + rise + rhythm,
          level + rise, level; plus a curve-level intercept in within designs)

- Starting values: the trajectory model's own cell values at its tau.
- Covariance: `vcov()` times N / (N - p) -- nlme's `adjustSigma`, which its own
  `summary()` and `anova(..., L = )` apply to ML fits. With it the Wald F here
  equals `anova(fit, L = )` exactly (checked in `test-traj-sat.R`).
- Denominator df: participants minus between-participant cells, for every
  test and interval. **In-house choice**: nlme assigns the within-participant
  df to every coefficient of a nonlinear model (597 on 48 participants x 14
  times), far too generous for comparisons between groups of participants; this
  is the between-within rule lme applies to between-participant terms.
- tau: Wald on log tau, interval back-transformed (geometric, asymmetric in
  hours); pairwise as a ratio. A: Wald. A/tau (the initial slope, the trend's
  rate at the first observation): delta method from the joint covariance of
  (A, log tau), gradient (1/tau, -A/tau) -- checked against a numerical
  derivative.
- "Is one tau enough?": the same model with one log tau, fitted by ML on the
  same random structure; likelihood-ratio chi-square on (cells - 1) df. Equals
  nlme's own `anova(fit_common, fit)`.
- Trend curves per cell: A (1 - exp(-t/tau)) with a delta-method pointwise band.
- Flags: |corr(A, log tau)| > 0.9 (A and tau not separately identified), and
  tau longer than the recording (A is an extrapolation).

**Two-stage**: each participant's A/tau (`init_slope`) from their own A_sat and
tau, compared like every other per-participant scalar (ANOVA, Welch pairwise);
group trend curves are the means of the participants' own A (1 - exp(-t/tau)).

**The factor in the change's unit line** is the trend column at the last
observation: T hours (linear), log(1 + T) (log), 1 - exp(-T/tau) (saturating:
the share of A reached by the last observation -- 0.8162 for T = 31 h and
tau = 18.3 h). It is now written out with T and tau.

Tests: `tests/testthat/test-traj-sat.R` (51 expectations: recovery, nlme
agreement, delta method, curves, flags, profile), `tests/traj_sat_app_test.R`
(the Run button under each tau option, per cell, flat profile, two-stage).

### tau: one value, or a range every estimate stays inside

Set with the saturating trend under Model Specification, and applied to every
fit that estimates tau:

- **Mixed, "Estimate one tau within a range"**: the profile grid is 24 points
  evenly spaced in log tau between the range's ends (it used to run from
  span/50 to 3 x span whatever the user wanted). A maximum at an end of the
  range is reported as a bound, not an estimate.
- **Mixed, per cell (nlme)**: tau is carried on a logistic scale,
  tau = lo + (hi - lo) * plogis(et), so the search cannot leave the range --
  the standard reparameterisation for bounding a parameter of a model without
  box constraints. Tests and intervals for tau go through log tau by the delta
  method (d log tau / d et = (hi - lo) dlogis(et) / tau); intervals are mapped
  back from the et scale and lie inside the range. A cell whose data put tau
  outside the range drives et to infinity and the fit fails; one unconstrained
  refit then names the cell and where its tau lies, as a diagnosis only.
- **Two-stage, "Estimate tau for each participant within a range"**: the range
  is the nonlinear least-squares bound on every participant's tau -- including
  the optimiser's last retry, which used to drop every bound for the saturating
  trend (the source of taus of hundreds of hours). Fits that end on a bound are
  flagged and counted, and "Include fits that hit a parameter bound" decides
  whether they enter the summaries. Previously tau was bounded by the user only
  behind "Enable Parameter Bounding" (which also bounds the intercept and the
  amplitudes), and otherwise ran from 0.5 h to five times the recording.
- **Two-stage, "Hold tau at a value"**: tau is a constant of the formula, A is a
  least-squares estimate; the trend-only null of the zero-amplitude F test holds
  the same tau (closed form), the parameter count excludes tau, and tau leaves
  the MANOVA, the component tests and the pairwise menu (a constant has no
  variance to test). The free-vs-fixed AIC check and the ridge diagnostics do
  not apply and say so.
- The nested-model (Delta-AICc) table fits its saturating candidates with the
  same range or held value, so its "reported" row is the fit reported.

Tests: `test-traj-sat.R` (the range bounds every per-cell tau and interval; the
constrained fit equals the free one when the optimum is inside; A/tau and the
curve band through the logistic map against numerical derivatives; a range that
excludes a cell's tau fails with the cell named), `traj_sat_app_test.R` (the
profile stays in the range; an end-of-range estimate is called a bound; an
impossible range is refused; two-stage held tau and ranged tau end to end).
## 7. Assumptions audit: nothing a run decides is left unstated

Every choice a cosinor run makes that is not a setting on screen -- a default,
a convention, a rule of the estimator -- is now listed for that run, read off
the fitted objects, by `dance_cosinor_assumptions()`
(`server/08j_helpers_assumptions.R`). The summary prints it as **Assumptions
and fixed settings** (and so does the summary download); the APA report puts
the same list in its Methods. In-house: only the wording; every rule it states
is the code's.

What the list covers:

- **Time**: where the times came from (clock times parsed at import, typed,
  a column, or -- the "column index" choice -- NOT read from the data but
  taken as equally spaced over exactly one period); the midnight rule; t = 0
  at the recording's first time point for every participant and cell; the
  period fixed, not estimated; cumulative harmonics.
- **Data**: raw observations; a missing one is left out, never imputed.
- **Trend**: log(1 + t) with t in hours (the +1 h is a convention); the MESOR
  window (first to last time point) and the change (the trend at the end).
- **Mixed**: the basis crossed with every design factor; REML via
  lme4/lmerTest, bobyqa; random effects normal with full covariance per term,
  the ladder rule (down only on failure or non-convergence, singular fits
  kept) and the rung used; **residuals normal, independent, one variance --
  no autocorrelation and no unequal variance is modelled**; the df method;
  marginal effects with equal weights; Wald-z per-cell intervals; joint
  amplitude/acrophase draws (6,000 per cell, 20,000 per contrast, seed
  20240601) with Bingham's rule; band types; tau's procedure (held value, or
  the profile's grid, structure, refinement, 1.92 cut, flat rule, bound rule,
  conditional inference); the per-cell nlme settings.
- **Two-stage**: the fitter (closed-form OLS; for the saturating trend
  nlsLM -> nls port -> four restarts, starts and the range), default bounds
  when bounding is on, the inclusion rule (non-converged out; on-a-bound
  in or out; the fewer-than-3 fallback), group summaries (arithmetic means;
  amplitude-weighted vector means; groups under 3 not summarised; unlabelled
  participants out of group analyses), the tests and their assumptions
  (estimates treated as exact; MANOVA Wilks; ANOVA equal variances;
  Bingham; Watson-Williams; Welch + Holm), rhythm detection at p < .05
  unadjusted, LOOCV, the bootstrap, the nested set and the free-vs-fixed
  check.

Behaviour the audit found wrong, and fixed:

1. **Clock times parsed at import were not unwrapped at midnight.** They are
   hours of the day in [0, 24), and the "shared clock times" route used them
   as they came, so a recording from 08:00 to 15:00 the next day put its second
   day BEFORE its first: every trend was fitted to scrambled times and t = 0 was
   midnight. All three clock routes now go through `dance_unwrap_clock()`
   (`server/03_helpers_clock.R`): a time earlier than the one before it is a
   midnight crossing and it, and every later time, moves 24 h on.
2. **The hand-entered (and column) times were unwrapped by one PERIOD**, not
   24 h -- the same thing only when P = 24 h. A clock does not wrap at noon
   because the period is 12 h.
3. **A participant missing the first observation had their own trend origin.**
   The two-stage fitters anchored the trend at each series' first NON-MISSING
   time, so that participant's intercept was the midline at a later hour than
   everyone else's, and group means averaged values from different times. Both
   fitters take `t_origin`; the run passes the recording's first time point, so
   every participant shares it (the mixed model is anchored there explicitly
   too, even if nobody was observed at it). The exported script passes it as
   well, and passes a held tau, which it did not.
4. **The nested-model table averaged each candidate over whoever it fitted**,
   while saying "same participants in every cell". `dance_model_set_common()`
   keeps the participants scorable under every candidate (a candidate fitted to
   nobody is named and dropped), and the summary and report state how many.
5. **Free tau against held tau compared two different samples** (free over the
   summarised fits, held over everyone with enough points), while saying
   "scorable in both". It is now a per-participant difference on the same
   observations and the same likelihood (`dance_gaussian_loglik`), averaged over
   the participants with both.
6. **The two-stage bootstrap was unseeded**: the same Run gave a different
   interval each time. It is seeded (20240601), the session's random stream is
   restored afterwards, and the summary states "percentile ... seed 20240601".
7. **A flat tau profile fell back to a hidden value.** The value it falls back to
   is the held value, now shown beside the range under the mixed approach in
   both modes, and it must lie inside the range (refused otherwise, as typed and
   at Run).
8. **Converged nonlinear fits were labelled "not_attempted"**: the status was
   assigned with `<<-` from inside `tryCatch`, which wrote it to the enclosing
   environment. The fit's own `converged` flag was always right.
9. **A time column of the wrong length silently became equally spaced times.**
   It is refused.
10. A group with no usable A or tau got a substituted A = 0 or tau = 1 h for its
    curve; it now has no trend curve (NA). The ridge diagnostics' correlation
    average no longer stops at the first 200 participants. The tau profile
    counts the grid points it left out. The long-recording note said values
    "will be wrapped using modulo"; nothing is.

Also: the duplicated tau block under Model Specification (an editing error in
the previous change, which created five input ids twice -- Shiny binds one and
the other does nothing) is gone, and a test renders the tab and fails on any
repeated id.

Tests: `test-cosinor-assumptions.R` (midnight unwrapping by 24 h on every route;
a missing first observation keeps the shared origin and recovers the true
intercept; converged fits say so and leak nothing; the nested set's common
sample; the assumption lists for mixed, flat-profile, failed-fit and two-stage
runs; the printed width), `test-harmonic-ui-structure.R` (every id on the tab
created once), `assumptions_app_test.R` (through Run: the clock route, the
shared origin, the refusals, the seeded bootstrap, the nested set's sample, the
summary and the report's list).

### The exported script reruns the cosinor analysis

The code export (`server/90_export.R`, section 10) used to refit each
participant's cosinor and stop. A mixed-effects run -- the default approach --
was not reproduced at all: the script printed per-participant estimates where
the app had reported a mixed model, and said nothing. The section also read its
data as `raw_data`, which the script never defined (section 2 creates
`data_matrix`), so for any cosinor fitted to raw observations -- every cosinor
since smoothing was removed from it -- the script stopped at that line. The
function-on-scalar section had the same hole.

Now the section reruns whichever approach ran, with the run's settings read off
the fitted model, by calling the functions the app calls -- moved out of the
module's server code so that there is one copy:

- **Mixed**: `dance_traj_run()` (server/08e) -- the long frame to the fitted
  model, tau held or profiled within the range (flat profile -> the held value),
  A and tau per cell on request -- and `dance_traj_tests()` (server/08f), every
  design test the comparison tab reports. The script prints the fit report, the
  fixed effects, the per-cell equations and intercept/MESOR/change, the design
  tests, and the per-cell A, tau, A/tau.
- **Two-stage**: `dance_ts_fit_all()`, `dance_ts_gate()`, `dance_ts_bootstrap()`,
  `dance_ts_model_set()` (server/08h): the per-participant fits, the inclusion
  rule, the seeded bootstrap and the nested set; then the comparison tab's
  tables and tests through the existing `dance_ts_*` functions.
- The run's assumptions head the section as comments (`dance_cosinor_assumptions`).
- Helpers are pulled in through codetools when it is installed and, when not,
  through every symbol in a function's body (a superset, narrowed to the app's
  own objects); app constants such as `DANCE_TRAJ_TRENDS` are written out too.
  Without the fallback a missing codetools emitted the requested function alone.

Not rerun by the script, and said so in it: the plots, and the two-stage ridge
diagnostics. Pairwise tables are one call each (`dance_traj_contrasts`,
`dance_ts_pairwise`), shown as comments.

Test: `tests/cosinor_export_test.R` presses Run, writes the script, runs its
header, data loading and cosinor section in a separate `Rscript --vanilla`
process with nothing of the app loaded, and compares: mixed (tau profiled, A
and tau per cell) -- tau, fixed effects, per-cell levels, all design tests' F,
df and p, per-cell A/tau/A-over-tau, the one-tau likelihood ratio; two-stage
(tau in a range, groups, bootstrap, nested set) -- the whole parameter table,
the inclusion audit, group summaries, component tests, bootstrap intervals and
AICc. Every number matches exactly.

## 8. Autocorrelated residuals in the mixed cosinor (continuous-time AR(1))

**The model.** Under Model Specification the mixed approach offers
"Residuals within a curve: Independent / Autocorrelated: continuous-time AR(1)".
With the AR(1), the residuals of one curve correlate as

    corr(e(t), e(t')) = phi^|t - t'|      (t in hours, 0 <= phi < 1)

-- nlme's `corCAR1` (Jones 1993; Pinheiro & Bates 2000, sec. 5.3), which takes
uneven spacing and missing observations as they are. Residuals of different
curves stay independent, with one variance throughout.

**How it is fitted: whitening, so the validated tests still apply.** The CAR(1)
correlation matrix has a bidiagonal inverse Cholesky factor W (Prais & Winsten
1954, for unequal gaps): along each curve in time order,
e*_1 = e_1 and e*_k = (e_k - r_k e_(k-1)) / sqrt(1 - r_k^2), r_k = phi^(t_k - t_(k-1)).
W R W' = I exactly, so for a given phi the CAR(1) model IS the ordinary mixed
model fitted to W y, W X and W Z with independent residuals (generalized least
squares), and the log-likelihood of the data -- ML or REML -- is the whitened
one plus log|det W|. The fit uses lme4's own modular steps (lFormula,
mkLmerDevfun, optimizeLmer, mkMerMod) on the whitened design, so the model
keeps its formula and coefficient names, and Kenward-Roger (pbkrtest) and
Satterthwaite (lmerTest, from the same whitened deviance function) apply
unchanged, conditional on phi. In-house: the whitener (`dance_car1_whitener`)
and the fit wrapper (`dance_traj_lmer_car1`), server/08e.

Checked against nlme (`tests/testthat/test-traj-car1.R`): at a given phi the
whitened fit has the fixed effects, standard errors, sigma and REML
log-likelihood of `nlme::lme(correlation = corCAR1)`; the phi profile finds
nlme's REML estimate; its likelihood ratio against independence is nlme's;
phi = 0 is the independent fit exactly. The same holds with a
within-participant factor, where the correlation runs within each curve
(participant x level) -- with a curve-level random effect, and with random
effects by participant only (nlme's `corCAR1(form = ~ t | subject/curve)`). nlme itself was not used for the
tests: its default ("containment") df ignore random slopes -- 200+ df for a
between-group rhythm contrast where Kenward-Roger and Satterthwaite give about
18 -- and emmeans' Satterthwaite for lme objects differentiates the covariance
by refitting on randomly perturbed data (unseeded), so its df change between
identical runs.

**How phi is estimated.** By profile likelihood on the random structure the
independent fit reached: 20 values from 0 to 0.98, evenly spaced in
-log(1 - phi), the maximum refined between its neighbours, the 95% interval
where the log-likelihood falls 1.92, and the likelihood-ratio test of phi = 0,
whose p is half the chi-square(1) tail because 0 is on the boundary (Self &
Liang 1987). REML when the fixed-effect design is fixed (tau held, or a
linear, log or no trend); maximum likelihood when tau is estimated too, since
REML likelihoods of different trend shapes are not comparable -- then tau and
phi are profiled in turn, each given the other, until phi moves by less than
0.005 (at most four rounds). Then phi's interval is the profile at tau-hat (tau
enters the mean and phi the covariance, and in a Gaussian model their expected
cross-information is zero, so to first order this is the full profile
interval; the reports say it is computed at tau-hat), and the likelihood ratio
compares with the independent model at ITS best tau -- the maximum of the
first tau profile, which is computed at phi = 0 (profiled again on the final
random structure if the rounds moved to another), not the independent model at
the tau the AR(1) chose, which would overstate the ratio. The model is then
refitted by REML at phi-hat; its tests treat phi (like tau) as known.

**What changes downstream.** R-squared and the residual plots are computed on
the data's scale (the whitened model's own `fitted()` is not); the residual
plots show the normalized (whitened) residuals; the information criteria use
the data's log-likelihood and count phi; the per-cell nonlinear model (nlme)
gets corCAR1 too, with phi estimated within it (nlme 3.1-164 crashes R when a
corCAR1 is held fixed in `nlme()`). The summary, the comparison tab's header,
the assumptions list, the APA Methods and Results and the exported script
state the residual model, phi-hat, its interval and its test. For an
independent-residual fit the summary reports how strongly consecutive residuals
of a curve correlate -- with no threshold attached; the AR(1) model's
likelihood ratio is the test.

**Calibration** (`tests/traj_car1_calibration.R`; two groups of 12, 16 uneven
times, one harmonic, random intercept and rhythm per participant, a few
missing points; circadian block test of Group, Kenward-Roger, alpha = .05):

| data (300 replicates each)        | independent-residual model | AR(1) model | phi-hat (true) | phi interval covers |
|-----------------------------------|----------------------------|-------------|----------------|---------------------|
| null, independent residuals       | 0.043                      | 0.040       | 0.104 (0)      | 0.987               |
| null, AR(1) residuals, phi = 0.6  | 0.050                      | 0.037       | 0.589 (0.6)    | 0.947               |
| 0.75 h phase shift, AR(1), power  | 0.687                      | 0.670       | 0.587 (0.6)    | 0.947               |

Read plainly: the AR(1) model holds its size (slightly conservative), estimates
phi without material bias, and its profile interval covers at the nominal rate
(above it when phi = 0, which is the boundary); power is essentially the same.
In this design -- a random rhythm per participant -- the independent-residual
model's between-group test held its size too: the random effects absorb most
of what independence gets wrong for a between-participant contrast. What the
AR(1) adds is a correct residual model, with phi estimated and tested; it was
not measured for designs whose random structure falls back to fewer terms,
for within-participant factors, or for other sampling grids.

**Fixed along the way.** (1) The per-cell nonlinear model grouped its CAR(1)
by participant on the rung whose random effects are by participant only; with
a within-participant factor that joined a participant's sessions into one
series (and nlme refused the repeated times). It is grouped by curve on every
rung now. (2) A fit restarted on an already-chosen random structure below the
first rung -- what the tau and phi profiles do -- stopped with an error while
listing the rungs it had skipped, so a tau profile on a reduced structure
failed outright. The skipped rungs are now recorded as "not refitted
(structure carried over)" in the fit's log. (3) phi grid points that do not
converge are left out of the profile and counted in the residual statement,
as for tau.

## 9. "No smoothing" is the data points joined by straight lines

The fPCA family (fPCA, time-warped PCA, fANOVA, clustering) needs an fd object
even when Data Preprocessing is set to "Raw data (no smoothing)", or when no
smoothing step was run. That object was a cubic B-spline with as many functions
as time points on fda's default, EVENLY spaced knots
(`create.bspline.basis(nbasis = n_time, norder = 4)`). It passed through every
observation, but its knots did not sit at the observation times, and between
them it swung far outside the data -- worst next to the first and last
observations. On 18 hourly VAS sleepiness ratings whose means run from 10 to 65,
the fPCA mean curve ran from -2 to 79 (reproduced with base R's B-splines in
`tests/testthat/test-fd-no-smoothing.R`), so it did not look like the mean of the
raw data.

Now (`dance_fd_no_smoothing`, server/04_helpers_fd.R): an order-2 B-spline basis
with a knot at every observation time. Each function is 1 at its own time point
and 0 at the others, so the coefficients ARE the observed values -- no fitting
step -- and each curve is the line plot of its data; the mean curve is the line
through the column means. A missing value, which an fd object cannot carry, is
read off the straight line between the curve's neighbouring observations (held at
the nearest observation before the first or after the last); it used to be the
participant's own mean. This applies to that representation only: the raw matrix
keeps its NAs. The APA report, the smoothing summary and the exported script say
so (the script previously built no curves at all for a no-smoothing run).

The same holds after time warping. The registered curves go back onto the
curves' own basis through one helper (`dance_refd`): straight lines through the
registered values when nothing was smoothed, the smoothing basis otherwise (as
the registration kernels already did). Time-warped PCA used to put them on
min(20, n - 2, max(4, n / 3)) cubic B-splines of its own -- 6 functions on 18
time points, a smoothing step nobody asked for -- and fANOVA and its pairwise
tests on min(20, n - 2), when the kernel returned no fd object.

## 10. How much warping, and is it needed?

There is no GCV for the amount of warping: every registration criterion
improves as the warps get more flexible, down to absorbing amplitude
differences ("pinching"); Marron, Ramsay, Sangalli & Srivastava (2015,
*Statistical Science*) leave the choice to the scientific question. Three
things can be estimated, and the time-warped PCA results now report them
(server/05b_helpers_warpstats.R).

**Is warping needed? -- the shape-invariant model** (`dance_shift_model`). One
common curve for everyone, moved up or down and shifted in time per
participant: y_ij = b0 + a_i + f(t_ij - s_i) + e_ij, a_i ~ N(0, sd_a^2),
s_i ~ N(0, sd_s^2), f a natural cubic spline (df set on screen, default 6;
it continues as a straight line beyond the data, where a shifted time can
land). Self-modelling regression (Lawton, Sylvestre & Maggio 1972; Kneip &
Gasser 1988; Brumback & Lindstrom 2004). Fitted by nlme, maximum likelihood,
on the RAW data in the data's time units; the model without shifts is linear
and fitted exactly by lme. Reported: sd_s with its 95% interval
(nlme::intervals), the likelihood ratio against no shifts (sd_s = 0 is on the
boundary: half the chi-square(1) tail), AIC and BIC, and each participant's
predicted shift (a random-effect prediction, shrunk towards zero).
Simulated (24 participants, 18 hourly points, a rising sigmoid, 200 data sets
each): with no shifts the test rejected 3.1% at .05 (fitted 193/200; the rest
failed with the spread at zero, which the message says); with sd_s = 0.75 h it
rejected 85%, mean estimate 0.71 h. Shift only -- no stretching of time --
and independent residuals.

**Amplitude and phase** (`dance_amp_phase`): Kneip & Ramsay's (2008, *JASA*)
decomposition, computed by fda::AmpPhaseDecomp from the curves before and after
registration and the warping functions (as straight lines through the warp
values; for a translation the identity, whose derivative -- all fda reads --
is the translation's). It needs warps that map the recording onto itself, so
it is reported for parametric and landmark warps and periodic shifts, and
refused for a shift of a non-periodic recording, where it is not a
decomposition (on test data it attributed -107% of the variation to phase).

**Each participant's warp, by group** (`dance_warp_estimates`,
`dance_warp_group_test`): the shift in the data's time units (positive = the
participant ran later than the reference), the parametric alpha with the mean
displacement of the time axis, or that displacement alone (landmark). Groups
compared by Welch's one-way ANOVA (stats::oneway.test) and Kruskal-Wallis
(stats::kruskal.test); periodic shifts are angles, compared by the
Watson-Williams test already used for acrophases.

**Fixed: the non-periodic linear shift registration.** Its lag came from ccf(),
which divides every lag's sum by the full series length and so pulls the
estimate towards zero -- strongly for curves that rise across the recording (a
sigmoid two steps late scored 0.96 at lag 0, 0.83 at its true lag) -- and the
shift was applied with the wrong sign, moving every curve away from the
reference: identical peaks at 0.4, 0.5 and 0.6 registered at 0.4, 0.6 and 0.8.
It also divided by n_time where a grid step is 1/(n_time - 1). Now: the
correlation of the curve with the reference over the times both cover, at each
lag up to a quarter of the recording, the best lag applied with the sign that
aligns, ties to the smallest lag. Time-warped PCA and fANOVA runs that used
the non-periodic linear shift before this fix registered the curves wrongly.
The periodic shift was correct.

**In the APA report and the exported script.** Methods: the shift model (as
above), the group comparison of the per-participant registration estimates,
and the amplitude/phase decomposition where it is defined. Results: SD(s) with
its 95% interval, the likelihood ratio, AIC/BIC; per-group means and SDs of
the estimates with Welch's F, Kruskal-Wallis H (or Watson-Williams), and the
decomposition. The script writes out the same functions and calls them with
the run's settings (the time axis as numbers); `tests/warp_report_export_test.R`
runs it in a clean `Rscript --vanilla` session and finds SD(s), its interval,
the test, the predicted shifts and every participant's warp identical to the
app's.

## 11. Several curves per participant (repeated cycles, sessions)

Found on a repeated-cycles file (13 participants, up to 7 cycles each, 85
curves of two-hourly sleepiness ratings with gaps):

- **The design readout identified participants by ROW**, not by the
  participant column read at import, so every curve was its own participant
  and a within-participant factor ("cycle") read as between. The fit itself
  already used the participant column; the readout now does too, and says when
  no identifier was read.
- **Two-stage stopped with "subscript out of bounds"** when a curve had no
  joint confidence region for a harmonic (a participant at the floor of the
  scale, 1 at every time point, has zero amplitude). The region is NULL, and
  `bingham[[h]] <- NULL` deletes the list slot instead of storing NULL, so the
  regions shifted and the second harmonic's was out of bounds. Stored with
  `bingham[h] <- list(...)` now.
- **The p-value formatter** coerces its input to one number and calls
  base::sub explicitly; the reported failure came through a sub() that was not
  base R's, presumably from a package attached in the user's session.
- **The data status** counts curves and participants ("85 curves x 14 time
  points; Participants: 13 -- repeated measures: each contributes 3 to 7
  curves").
- **The raw-data mean curves** (data plot, group previews) are means of the
  OBSERVED values at each time point; without na.rm a single missing value
  removed that time point from the mean.

A within-participant factor needs the mixed-effects approach: two-stage
compares the per-curve estimates as independent groups and says so.
`tests/repeated_curves_app_test.R` covers the case.

## 12. Two-stage estimated trajectories: the trend views

With a trend in the model, the two-stage "Estimated trajectories" offered only
the full trajectory and the rhythm with the trend removed, so differences in
the trend between groups could not be drawn. `dance_ts_group_curves()` now
takes the component the mixed approach offers -- full, the non-periodic change
from the first observation (trend(t) - trend(t0), zero at the start), intercept
+ harmonics, harmonics alone -- and computes each group's line from its
vector-mean coefficients and the band from the participants' own curves, as
before. Overlaying the raw data under two-stage failed with "object 'tv' not
found" (a variable from another function); fixed.
`tests/ts_trajectory_views_test.R`.

## 13. A package in the user's session that masks a base function

The cosinor summary failed with "subscript out of bounds" inside sub(). The
traceback ran sub -> .get_parts -> isFALSE, which is not base R's sub: a package
attached in the user's R session before the app started exported its own `sub`
(rsat is one CRAN package that does), and every sub() in the app reached it.
app.R now records the search path when it starts, and binds base R's version
of every function such a package masks in the environment the server functions
look in first, and says so at startup ("DANCE: using base R's sub() (masked by
...) inside the app"). The user's workspace is not changed. DANCE's own packages
and their dependencies are exempt: Matrix, attached with lme4, extends t(),
solve(), as.matrix() and others for sparse matrices, and the mixed model relies
on those. `tests/base_masking_app_test.R` loads app.R as runApp does with a
masking `sub` attached; with the old app.R it fails exactly as reported.

The same can happen to DANCE's own packages. "Group Comparison on Component
Scores" failed with "is.character(txt) is not TRUE": that is the first line of
jsonlite::validate(). When shiny is already attached in the session and jsonlite
is attached after it, library(shiny) in app.R does not move shiny, so jsonlite's
validate() sits in front of shiny's and every validate(need(...)) in the app
fails. app.R now also binds, for every function a package outside DANCE's set
masks from one of DANCE's attached packages, the version a fresh session would
find (the first of DANCE's packages on the search path that exports it), and
says so at startup ("DANCE: using shiny::validate() (masked by jsonlite) inside
the app"). The same test attaches shiny and then jsonlite first.

Also: a column named ..._24 is read as midnight when the recording carries on
past it (22, 24, 2, ...), so CORT_10 ... CORT_24, CORT_2 ... CORT_8,
CORT_10_SD ... CORT_20_SD give 36 h of clock time; a series that merely ends
at 24 (T1 ... T24) is still not taken for clock time.

## 14. Two-stage estimated trajectories: the mixed model's views and bands

The two-stage panel now offers what the mixed one does, where a two-stage
counterpart exists:

- **Each harmonic on its own** (with more than one harmonic): harmonic k's
  cosine and sine coefficients alone, centred on zero.
- **Simultaneous bands** for the group curves and for the difference between
  two groups. When the view is a fixed linear combination x(t)'b of the
  coefficients (no trend, a linear or log trend; any view without the trend),
  each participant's b_i is one draw, and Hotelling's T^2 gives a band that
  holds for the whole curve at once: mean curve +/- c SE(t), c^2 = d (n - 1) /
  (n - d) F(d, n - d), d the coefficients the view uses, SE(t) the SD of the
  participants' curves at t over sqrt(n) (Scheffe's method with Hotelling's
  T^2). The difference of two groups: the two-sample version with pooled
  covariance, c^2 = d (N - 2) / (N - d - 1) F(d, N - d - 1). Simulated (2000
  data sets each, d = 6): the whole true curve was covered in 1.000, 1.000 and
  0.995 of them for n = 10, 20, 40 -- conservative, as Scheffe bands are --
  against 0.66 to 0.75 for the pointwise band. With small groups the band is
  wide (n = 12, d = 6: c = 6.86). With a saturating trend every participant
  has their own tau, the curves are not linear in shared coefficients, and the
  band stays pointwise; the note under the plot says so.
- The difference plot takes the same pointwise/simultaneous choice as under
  the mixed approach (`dance_ts_diff_curve`).

`tests/ts_trajectory_views_test.R`.

## 15. Mixed cosinor pairwise: the trend coefficient beside the change

The pairwise menu offered the non-periodic change over the recording but not
the trend coefficient that the model summary reports (the saturating rise A, the
slope, or the log-trend beta). The two are one hypothesis: the change from the
first to the last observation is the coefficient times a known constant,

    change = A * (1 - exp(-T / tau))      (saturating; T = last - first observation)
    change = beta * T                     (linear)
    change = beta * log(1 + T)            (log)

so the F, the df and the p-value are identical, while the estimates differ by
that constant. The coefficient is now offered as its own entry ("Saturating rise
A (the asymptote above the intercept)", "Slope beta (per hour)", "Log-trend
beta"), immediately before the change, and its unit line states the factor that
turns it into the change. A is the rise the trend approaches as t grows. It is
not the rise at t = tau, which is A * (1 - 1/e), about 63% of A. When the
per-cell nonlinear model (a tau per cell) has been fitted, its A replaces this
entry, so the menu never shows two different numbers called A.


## 16. An empty participant identifier

A file whose identifier column has an empty cell (one row in a real file: the
cell between subjects 22 and 24) passed NA on as that participant's
identifier. The participant picker showed the row as "Row 22", and the mixed
cosinor dropped all of its observations, because lme4's model frame removes
rows whose grouping factor is NA (on a 24-participant test file: 23 groups,
184 of 192 observations). The two-stage approach and the other tabs kept it.
Each row of a wide file is one curve, so an unnamed row is one participant
of its own. The import now names it "Row N (no ID)" (N is its row of the data,
counting from the first row under the header), says so in a notification and
in the data status, and every analysis sees the same participants.
`tests/missing_id_app_test.R` fails on the old import (23 participants in the
mixed fit) and passes on the new one.

## 17. The cosinor's pairwise tables in the exports

Every quantity the "Compare" menu offers, mixed or two-stage, is written to
one CSV (under the pairwise panel, and on the Export tab). It covers the whole
design, with the multiplicity adjustment chosen on the tab, and its family is
the pairs of one quantity, as on screen. Columns: the quantity and its unit,
the pair, the estimate with its interval, the statistic and its df, p and
adjusted p, and Cohen's d for the two-stage Welch tests. Joint block tests
have no single estimate, and their rows say so (joint_test = TRUE). The
exported R script computes the same tables (traj_pairwise or ts_pairwise)
instead of naming the function in a comment. `tests/cosinor_export_test.R`
compares the script's tables with the app's.

## 18. Header cells formatted as times

Excel stores a time as a fraction of a day (08:00 = 1/3). A header row of
8:00, 10:00, ..., 0:00, ..., 16:00 formatted as times reached the app through
readxl as "0.33333333333333331...5", "0.75", "0", "8.3333333333333329E-2",
..., with "...N" added to the times that recur, and the hour parser read 0.333
h. At import such names are converted back to the clock time they display
(minutes = fraction x 1440, rounded to the minute), and a recurring time is
numbered ("08:00 #2") so that the names stay unique. The conversion is applied
only when it is clearly what happened: at least two names that are numbers in
[0, 1) falling on a whole minute, at least one of them with the full binary
precision Excel writes. Typed decimal hours such as 0.25 or 0.5 are left as
they are. The import says which columns it converted.

The hour parser also reads "HH:MM" with its minutes. It used to take the first
number in the name, so 8:30 was 8 and MT_21:30_SD was 21.

## 19. Half-hourly columns, a gap, and where the plots put the hours

A protocol with one evening every 30 min (MT_20h00 ... MT_00h00), then the
next evening into a night of sleep deprivation (... MT_07h00, MT_16h00): 25
columns over 44 h, with a 20 h gap.

* The import read every column correctly (20, 20.5, ..., 0, 20, ..., 7, 16).
* The cosinor tab's "Specify times manually" suggestion ran its own column-name
  patterns. For MT_20h30 the "trailing number" rule took the minutes and
  suggested 0, 30, 0, 30, ... It now suggests the import's clock times, so the
  manual route and "Use shared clock times" start from the same numbers.
* Manual times are unwrapped at midnight: a time earlier than the one before it
  is the next day. "00:00 then 20:00" is the same day (+20 h), which is right
  for this protocol. A list typed out of order (23.5 before 23) adds a whole
  day at that point, and the fit then spans three days.
* The plots drew the curves at each column's position on the smoothing's own
  axis (the column index, unless real clock time was chosen for smoothing), but
  put the hour labels at REAL-time positions. On unevenly spaced columns the
  labels did not sit under the data. The labels and the hover times now use the
  positions the curves are drawn at (dance_column_positions()). The hover text
  also interpolates elapsed time, so a point between 23:30 and 00:00 reads 23:45
  and not 11:45.

## 20. Pairwise comparisons in the APA report

Both cosinor approaches report every pairwise table under "Pairwise
comparisons", built from the same function as the CSV export (section 17).
A one-df contrast is reported as t (z when the reference is normal, as for the
delta-method amplitude and acrophase contrasts of the mixed model), signed like
the difference. A mixed F(1, df) is reported as t(df) = sign x sqrt(F). Joint
block tests, and the Watson-Williams acrophase test of the two-stage approach,
are reported as F.

## 21. Three model structures, three scientific questions

### Model structures
A selector on the cosinor tab chooses one of three analyses:

* **rhythmic component only**: harmonics, no trend;
* **non-periodic trend only**: a linear, logarithmic or saturating trend, no
  harmonics;
* **both**.

The choice sets the trend and the harmonic count (0 for trend-only). Both
estimators, every output, the report and the exported script follow it. A
trend-only model has no amplitude, acrophase, rhythm test, polar plot or
rhythm-detection count; each place that would show one says why not.

### The outcomes, defined once (server/08_helpers_cosinor.R)
Time is in hours since the first observation, t in [0, T]. The trend is c·s(t)
with s(0) = 0:

| trend   | s(t)              | s(T)              | mean of s on [0, T]           | s'(0)  |
|---------|-------------------|-------------------|-------------------------------|--------|
| linear  | t                 | T                 | T / 2                         | 1      |
| log     | log(1 + t)        | log(1 + T)        | ((1 + T) log(1 + T) − T) / T  | 1      |
| exp_sat | 1 − e^(−t/τ)      | 1 − e^(−T/τ)      | 1 − (τ/T)(1 − e^(−T/τ))       | 1/τ    |

* total trend change ΔT = c·s(T)
* average rate = ΔT / T
* initial rate = c·s'(0)
* mean midline (MESOR) = b0 + c·mean(s)
* time-averaged fitted response = MESOR + Σ_h [β_h sin(ω_h T) + γ_h (1 − cos(ω_h T))] / (ω_h T)

None of the three trends has a singular derivative at 0: the log trend is
log(1 + t). For a straight line the initial and average rates are one number,
so only the average rate is offered.

These quantities are:

* linear functions of the mixed model's fixed effects (exact L-matrix tests);
* computed from each participant's own coefficients in the two-stage fits.

Checked against numerical integration and differentiation for every trend, at
every harmonic count, in both approaches (tests/testthat/test-outcomes-calculus.R,
agreement to about 1e-8).

### With τ estimated
* **Two-stage.** Each participant's ΔT, average rate and A/τ use their own A
  and τ.
* **Mixed, one shared τ (profiled).** The quantities are computed at the
  profiled τ, and the summary says the inference is conditional on it.
* **Mixed, τ per cell (nlme).** ΔT and the average rate carry τ's uncertainty
  through the delta method:
  ∂ΔT/∂A = 1 − e^(−T/τ) and ∂ΔT/∂x = −A (T/τ²) e^(−T/τ) dτ/dx, where x is
  τ's carrying scale.

### The comparison, by question (server/08k_helpers_outcomes.R)
Above the three domains sits the **omnibus** test of the complete fitted
trajectory (mixed: the effect tests; two-stage: the MANOVA on every
coefficient). Then three domains:

1. **Overall response level.** Time-averaged fitted response (primary), mean
   midline (MESOR), initial midline.
2. **Circadian-like rhythmicity.** The joint harmonic-vector test (primary;
   neither an amplitude nor an acrophase test), then amplitude and acrophase per
   harmonic.
3. **Homeostatic-like accumulation.** Total trend change (primary; one test
   with the average rate), initial rate, then the model-specific parameters
   (the trend coefficient, A, τ, t50 = τ ln 2).

Group intervals: Wald for the mixed model; t intervals of the participant means
for two-stage. Two-stage amplitude and acrophase intervals are Bingham et al.'s
(1982) population-mean region: dance_bingham_ci() with V = S(k − 1)/(k(k − 2))
on F(2, k − 2). No test was substituted for an existing one.

### Graphs
The trajectory plot's full, rhythm and trend views carry:

* the time-averaged fitted response line;
* the total change and the average-rate line;
* the initial tangent;
* the asymptote and the 63.2% point at τ;
* each harmonic's acrophase at its amplitude.

All are read off the same outcome table. tests/structures_app_test.R checks the
drawn values against it.

### Model comparison
The candidate space is intercept only, each trend alone, the rhythm alone and
each trend with each cumulative harmonic set, and it is printed in full.

* **Two-stage:** per-participant AICc, averaged over the participants fitted by
  every candidate.
* **Mixed (new):** each candidate fitted by ML (REML likelihoods are not
  comparable across fixed effects), with its own random-effects rung.
* **Parameter counts:** a profiled or per-participant τ counts as a parameter;
  a held τ does not.
* **Pairs that are not nested** (log against saturating) are compared by AICc,
  not by likelihood-ratio tests.

### Is τ identified?
Each approach gets a checklist:

* convergence;
* estimates at a bound;
* relative SE of τ;
* corr(A, τ);
* τ against the window;
* free versus held τ;
* the number of valid estimates.

The verdict is "identified", "weakly identified" or "not identified". Anything
short of identified recommends a fixed-τ sensitivity analysis. Hierarchical τ
(an advanced option) adds a participant-level random effect on τ to the
per-cell nlme model and reports its SD in hours.

### Compatibility mapping
Nothing exported was redefined: `mesor` (the intercept), `mesor_adj` (the
MESOR), `trend_change` and `init_slope` keep their meanings and columns. New
names and changed labels:

| where | before | now |
|---|---|---|
| labels | Intercept (midline at the first observation) | Initial midline (t = 0) |
| labels | MESOR (midline averaged over the window) | Mean midline level (MESOR) |
| labels | Non-periodic change over the recording | Total trend change (start to end) |
| per-participant columns | – | `fitted_mean`, `avg_rate`, `rate0`, `t50` (added) |
| pairwise keys (both) | – | `fitted_mean`, `avg_rate`, `rate0` (added) |
| pairwise keys, two-stage | `init_slope` | `rate0` (A/τ, the same number); `t50` added |
| pairwise keys, mixed per-cell | `A_sat`, `tau`, `init_slope` | `sat_change`, `sat_avg_rate`, `init_slope`, `A_sat`, `tau` |
| model labels | – | `<trend> only`, `intercept only` (labels with harmonics unchanged) |
| model comparison | trends × H1..HK (two-stage only) | + intercept only and trends alone; mixed by ML too; a `structure` column |
| session files | trend "none" = rhythm only | restored as structure "rhythm"; any other trend as "combined" |

## 22. The mixed model comparison follows the reported fit's residuals; the table's rows link to the graph

**Residuals in the comparison.** Before, `dance_traj_model_set()` fitted every
candidate with independent residuals, even when the reported fit used
continuous-time AR(1). It said so in its output, but the result was a comparison
of models other than the one being reported. Now, when the reported fit has
φ > 0, every candidate is fitted by ML at that same φ (the whitened fit of §8).
- **Log-likelihood:** the data's log-likelihood, that is the whitened fit's
  log-likelihood + log|det W| (`dance_car1_logdet`). This makes it comparable
  across candidates, which share W.
- **Parameter count:** φ counts as one parameter in every candidate.
- **What φ is:** φ is held at the reported value, not re-estimated per
  candidate. The note printed with the table (summary, APA report, script)
  says so.
- **Test:** `tests/testthat/test-model-comparison.R` checks k(φ) = k(0) + 1
  and checks the log-likelihood against a direct whitened refit.

**Rows linked to the graph.** Clicking a row of the "By question" panel does
three things:
- it selects that quantity in the Compare menu, when the menu offers it;
- it moves the graph to the view the quantity's annotation is drawn in (full
  fitted response, trend component, or the harmonic);
- it switches that annotation on.

Choosing another quantity in the menu by hand ends the row's selection.
- **Test:** `tests/structures_app_test.R`, all 16 structure × approach runs.

## 23. Dim-light melatonin onset (DLMO)

The new tab "Melatonin Onset (DLMO)" works on whatever data are selected at
Data Import: each row is a recording, each column a sample. The code is in
`server/08l_helpers_dlmo.R` (pure functions) and `server/76_dlmo.R` (the tab).

**Packages checked first.**
- **CRAN:** no package estimates DLMO.
- **GitHub:** **dlmoR** (Thalji & Spitschan 2026, *J Biol Rhythms* 41:301;
  MIT licence; tscnlab/dlmoR) implements the hockey-stick method of Danilenko
  et al. (2014, *Chronobiol Int* 31:349), and only that method.
- **How the hockey-stick runs:** it is called through dlmoR and is not
  reimplemented. Without dlmoR the method is refused with the install command,
  and the threshold methods still run.
- **Startup:** dlmoR is declared as an optional package to install from
  GitHub. It is never attached, because it exports generic names.
- **dlmoR's own default threshold** is 2.3 pg/mL. The app passes the user's
  threshold instead.

**In-house calculations** (no package does them; each is a textbook rule):

| Step | Rule | Source |
|---|---|---|
| Fixed threshold | Onset is the first sample *strictly above* the threshold whose next `confirm` samples (default 1) are above it too. The sample before it must be at or below the threshold. A rise that falls back sooner is skipped and counted. | Lewy & Sack 1989; Benloucif et al. 2008 |
| Interpolation | Linear, between the two **observed** samples that bracket the crossing. A missing sample is never filled in; the line spans its observed neighbours. If they are more than `max_gap` hours apart (default 2), the onset is not used ("gap"). | Benloucif et al. 2008 |
| Variable threshold | Mean + k·SD (default k = 2) of the first n **observed** samples of each profile (default n = 3). With n ≤ 5 these samples cannot themselves exceed it. | Voultsios, Kennaway & Dawson 1997 |
| LLOQ | Optional. A value below the LLOQ is set to LLOQ/2, the LLOQ, or 0 (the user chooses), and counted. Without an LLOQ, values are used as entered. | stated in the output |
| Phase angle | Bedtime − DLMO, in hours; positive when the onset precedes bedtime. The bedtime used is the one within 12 h of the onset. Computed **only** when a habitual bedtime is given (one time for everyone, or a column); otherwise there is no phase angle anywhere. | |

**Profiles and nights.**
- A gap longer than `gap_h` (default 6 h) between scheduled sample times
  starts a new profile. That makes "two nights in one row" two profiles,
  "Night 1" and "Night 2".
- With one row per night, a column names the night.
- Each profile is moved by whole days onto one evening clock: its first sample
  falls in [anchor, anchor + 24), where the anchor is the earliest start hour.
  So 01:00 is 25, after 23:00, and both nights are read on the same clock.
- Column names that read as the same time are refused, and the user is told
  to enter the times by hand: the hour comes from a name's first number, so
  "Mel_N2_18" reads as 2.

**Censoring: nothing gets a time it does not have.**
- **Before the window:** above the threshold from the first sample.
- **No onset:** never a sustained rise.
- **Unconfirmed:** a crossing too close to the last sample to confirm.
- **Gap:** the bracketing samples are too far apart.
- **Too few baseline samples** for the variable threshold.

Every one of these is reported and counted per cell, and none enters a mean
or a test. If being undetermined depends on the condition (an onset pushed
past the window after sleep deprivation, for instance), the comparison is
biased. The counts are shown so a reader can judge that.

**Comparisons** (`dance_dlmo_compare`), on determined onsets only:

| Factors | Test |
|---|---|
| Nights only, two nights | Paired t test on participants with an onset on both nights; the unpaired participants are counted. |
| Nights with groups, or more than two nights | `lmerTest::lmer(dlmo ~ group * night + (1 \| id))`: Satterthwaite F tests, and the night difference within each group as a model contrast (`lmerTest::contest`). |
| Groups only | Welch's t test (two groups) or Welch's one-way ANOVA. Refused if participants contribute several profiles. |

The phase angle is compared the same way.

When more than one method is run, their agreement is reported: the mean, SD
and largest absolute difference over the profiles where both determined an
onset.

**A property of the method, not a defect.** When the rise starts between two
sparse samples, the straight line between them crosses the threshold early.
For example, hourly samples of a rise starting at 20:30 put the 3 pg/mL
crossing at 20:40 instead of 20:50. This is why the onset depends on sampling
density (Molina & Burgess 2011). The test suite checks it explicitly.

**Tests.**
- `tests/testthat/test-dlmo.R`:
  - exact crossings, spikes, both censorings, unconfirmed, gaps;
  - variable threshold, LLOQ rules, night splitting, crossings after midnight;
  - phase angle and bedtime parsing, and dlmoR result conversion;
  - the hockey-stick refused without dlmoR, and checked against a known onset
    when dlmoR is installed;
  - recovery on simulated profiles, the paired and mixed comparisons with a
    known night × group effect, and method agreement.
- `tests/dlmo_app_test.R`:
  - both layouts through the app, with every output rendered;
  - the CSV, the APA report, and the exported script reproducing every onset
    and phase angle;
  - repeated times refused, and times entered by hand.

### 23a. The hockey-stick, tested with dlmoR installed (dlmoR 0.2.0)

- **A bug found and fixed.** dlmoR returns the onset as an `hms`, which is
  always held in seconds and ignores `as.numeric(units = "hours")`. The app
  was showing seconds modulo 24 as the onset, e.g. 11:00 where dlmoR said
  20:20. The time is now read from its stored value and unit.
- **A related fix.** A time of day is always placed on the day that falls
  within the sampled window, so a 07:00 onset in an 18:00–07:00 profile is
  the next morning.
- **Agreement with dlmoR.** Every profile dlmoR ships was run two ways: dlmoR
  called directly on the original date-times, and the app's route (hours on a
  reference date, read back). Fine search, threshold 5 pg/mL, interval limit
  2 h.

| profile | samples | dlmoR | app | fine-search time |
|---|---|---|---|---|
| sample | 19 | 20:20 | 20:20 | 150 s |
| b | 12 | 00:42 | 00:42 | 398 s |
| c | 12 | 20:07 | 20:07 | 246 s |
| d | 11 | error "no dynamic part" | status "failed", dlmoR's message as the note | – |
| e | 12 | 21:35 | 21:35 | 565 s |
| f | 9 | 20:06 | 20:06 | 70 s |
| g | 12 | 21:16 | 21:16 | 508 s |
| h | 12 | 19:42 | 19:42 | 187 s |
| i | 12 | 21:00 | 21:00 | 459 s |

  Every determined onset is identical to the minute, including b, whose onset
  falls after midnight.
- **Run time.** The fine search takes 1–10 minutes per profile, so the run
  shows progress per profile. The fine search can be switched off: coarse
  only, a few seconds per profile. The settings text states which search was
  used.
- **The app test** (`tests/dlmo_app_test.R`) runs the hockey-stick through the
  app when dlmoR is installed. It checks that the app's onsets equal dlmoR's,
  that they are within 1 h of the generating onsets, and that the report cites
  Danilenko and dlmoR. It also checks that the exported script reproduces the
  onsets.
- **What the hockey-stick estimates** on synthetic convex rises (a power 1.6
  rise starting at 20:30, 21:00 and 21:30, half-hourly samples, coarse
  search): 21:27, 21:51 and 21:57. These are 0.5–1 h after the start of the
  generating rise. The hockey-stick onset is where its fitted baseline and
  ascending segments meet, not where the curve first departs from baseline,
  so on a slowly accelerating rise it is later than the rise's start. It is a
  different definition from a threshold crossing, and the two are not
  interchangeable; the report says so.

## 24. The log scale, and intervals in which nothing was sampled (cosinor)

Raised on cortisol: a two-harmonic fit with no samples between 18:00 and 06:00
dipped to −1.7 overnight, which is impossible. Constraining the curve to the
observed range was considered and rejected:
- The true nadir lies inside the gap and is lower than any daytime sample, so
  that constraint would force a wrong nadir and bias the amplitude and the
  acrophase.
- An estimate held on a constraint has no valid standard error or test.
- There is no validated cosinor method with a shape constraint.

Two established remedies were added instead. Both are in
`server/08m_helpers_scale.R`.

**The log scale (the user's choice, never applied silently).** This is
"Scale of the analysis" in the settings, with a constant c (default 0).
- **What is fitted:** the model is fitted to ln(y + c). Every fit, test,
  bound and table is on that scale.
- **Non-positive values:** if any observed y + c ≤ 0, the run is refused
  with the count.
- **Read back on the original scale:** exp(x) − c.
  - The fitted curve cannot fall below −c.
  - A level (intercept, MESOR, time-averaged response) is a geometric mean
    when c = 0.
  - An amplitude A is the peak-to-trough ratio e^(2A) of its harmonic.
  - A trend change is a ratio e^(change), and a group difference in a level
    is a log ratio.
- **Intervals:** read back the same way. The transformation is monotone, so
  this is exact for the limits.
- **In the outcome tables:** the readings go in each row's note. The app
  panel, its CSV, the APA report (a Note column) and the exported script all
  carry them.
- **In the plots:** the fitted-curves tab and the full trajectory view show
  the original scale by default, with a checkbox to stay on the log scale.
  The finished figure is transformed as a whole: curves, bands, raw points
  (drawn as observed), reference lines, annotations and the hover value.
  Component views (harmonics, trend) are log ratios; they stay on the log
  scale and say so.
- **Admissible range:** checked on the original scale.
- **Parameter bounds:** these apply on the log scale, and the bounds hint
  shows log-scale values.
- **Stated** in "What this run assumed" (summary, report, script).

**Unsampled intervals.**
- **Definition:** a gap between consecutive sampling times longer than the
  user's threshold, or by default twice the median interval. Only times at
  which someone was observed count.
- **On the plots:** each gap is shaded and labelled "no samples", and each
  fitted curve is dashed across it. Raw-data lines are left solid.
- **Acrophases:** an acrophase is flagged when every occurrence within the
  window (every period/h for harmonic h) falls in a gap.
- **The fitted curve's maximum and minimum** in the summary are flagged when
  they lie in a gap.
- **Stated** in "What this run assumed".

**Three existing errors found on the way, and fixed.**
1. **Bingham amplitude interval.** When the joint (cos, sin) region contains
   the origin, it contains amplitude 0, so the lower bound is 0. The code took
   the smallest radius on the region's boundary instead. That point lies
   beyond the origin and gave a lower bound above the estimate itself, for
   example an estimate of 0.13 with an interval of [0.20, 0.59]. The
   acrophase was already, correctly, reported as not identified in that case.
   - Test: `tests/testthat/test-cosinor-audit.R`.
2. **Mixed summary extremes.** `dance_traj_curve_peaks()` returns
   `list(ok, table, ...)`, but the mixed summary read the list as the table.
   So it never printed the fitted curve's maximum and minimum. It now reads
   `$table`.
3. **Clock origin in the exported script.** The script's models carried no
   clock origin: the two-stage `ts_model` lacked `origin_shift`, and
   `traj_fit` lacked `clock_origin`. So every acrophase the script converted
   to clock time was off by the recording's start time (8 h for a recording
   starting at 08:00). The estimates themselves were right. Both models now
   carry it.
   - The round-trip test compares the clock times as well as the numbers.

**Tests.** `tests/scale_gap_app_test.R` uses cortisol-like data with an
18:00–06:00 gap and two groups, A peaking in the morning and B at midnight. It
runs both approaches on both scales and checks:
- the original-scale dip below zero, and that it is flagged;
- the log-scale fit on ln(y), with the read-back curve positive, raw data as
  observed, and the hover value matching the drawn value;
- the log-scale view on request;
- geometric means and peak/trough ratios that equal exp(x) and exp(2A);
- the gap shaded and the curve dashed only across it;
- B's acrophase flagged and A's not;
- the report;
- the exported script reproducing every estimate, note and clock time;
- a zero refused with c = 0, and accepted with c = 0.5.

## 25. Functional ANOVA: a curve without a group; the criterion the plots use

**The defect (found on real data).** A curve whose group value was missing
entered the between-groups functional ANOVA, with these consequences:
- **Observed F:** it left the curve out of the group sums, but counted it in
  the degrees of freedom (df within 28 instead of 27).
- **Permutations:** the curve's missing label was shuffled with the rest, so
  every permuted statistic came out undefined. As a result:
  - every pointwise p-value was missing: the p-value plot was empty and had no
    legend entries;
  - no region could be significant;
  - the L2 p was computed as (1 + 0)/(1 + 0) = 1, which reads as "no
    difference" when no comparison was made.
- **Reproduced:** on simulated data this gives exactly that screen, while the
  same data without the curve give L2 p = 0.005.

**The fixes.**
1. **Drop and count.** Both kernels (`perform_functional_anova`,
   `perform_rm_fanova`) drop a curve with a missing or blank group, subject
   or condition value before anything is computed, and count it
   (`dance_fanova_complete`). The readout, a notification and the APA report
   state how many. The result is identical to running without that curve
   (test).
2. **No permutation, no p-value.** `dance_perm_p()` returns NA when there are
   no valid permutations (B = 0), so the readout says "not computed". It
   never shows p = 1 in that case.
3. **Which criterion "significant" uses.**
   - The red line on the F plot is the uncorrected parametric critical F at
     one time point, and is now labelled so.
   - The stretches that are significant under the test's own criterion are
     shaded on the F plot, as on the mean-function plot: pointwise
     permutation p-values, FDR-adjusted (Benjamini-Hochberg) across the 100
     evaluation points.
   - The readout states that criterion.
   - An empty p-value plot says why.
4. **The exported script.**
   - **Labels:** it now passes the fANOVA tab's own grouping variable, with NA
     where a curve has no group or is outside the compared groups. It used to
     pass the first grouping variable selected at import, whatever the tab
     used, and wrote a missing group as the text 'NA', i.e. a group called
     "NA".
   - **Curves:** with no smoothing step, the app's analyses use the data
     points joined by straight lines (`dance_ensure_fd_obj`). The script
     fitted a 20-basis GCV smooth instead, so its functional analyses were on
     different curves. It now builds the app's representation.

**Tests.**
- `tests/testthat/test-fanova-missing-group.R`: the kernel and
  `dance_perm_p`.
- `tests/fanova_missing_group_app_test.R`: through the app, the readout, both
  plots, the report, and the exported script reproducing the same F.

## 26. DLMO: each group's mean or median onset on the profile plot

The profile plot's "Show" menu now offers all participants, each group (when
a group variable is defined) and each participant.
- **Summary lines:** with all participants or one group, a thick vertical line
  per method (solid fixed, dashed variable, dotted hockey-stick) marks each
  cell's mean or median onset (the user's choice, or none).
- **What a cell is:** a group × night combination for all participants; a
  night within one group.
- **What enters the summary:** determined onsets only. The values are the
  Comparisons tab's descriptives (`dance_dlmo_onset_summary`), and a table
  under the plot gives each value with "k of n profiles".
- This is display only: no number in the CSV, report or script changes.
- **Test:** `tests/dlmo_app_test.R`.

## 27. A streamlined comparison view; no "MESOR"; the equation behind each number

**The streamlined view is the default** of tab 6, Group / Condition
Comparison. It reads the same outcome table as the full output, so no number
in it is computed anywhere else:
- **Header:** the approach, the harmonics, the trend (with τ fixed or
  estimated) and, when chosen, the log scale.
- **Primary analysis:** the omnibus test of the complete fitted trajectory.
- **Secondary comparisons,** one row per question, each with its group
  estimates [95% CI] and its unadjusted test:
  - the overall response level (the time-averaged fitted response, or the
    time-averaged midline level for a trend-only model);
  - the non-periodic trend (the change over the window);
  - each harmonic's joint (cos, sin) test.
- **Rhythm descriptors (exploratory):** amplitude and component acrophase
  per harmonic. "Not estimable" is shown when the joint region contains the
  origin.
- **Quality checks:** only those that apply to the run. For example:
  - with τ fixed, the rate and asymptote tests repeat the change test;
  - the intervals are per group, not for differences;
  - the denominator degrees of freedom differ between tests;
  - bounds were applied;
  - the log scale;
  - unsampled intervals;
  - τ not identified.
- **Rows** select their quantity on the graph, as in "By question".

Everything else is kept, unchanged, under **Detailed statistics**: the
omnibus block, "By question", the component tests, the level, saturating-
trend and difference plots, and the pairwise tables. The estimated-
trajectories plot stays in view between the two.

**"MESOR" is gone from everything a reader sees:** the interface, tables,
plots, summary, downloads (column names and keys are written as
`midline_level`), the report, and the comments and printed labels of the
exported script.
- **The name used:** the time-averaged midline level, which is what the app
  computes. The midline (intercept + trend, harmonics excluded) is averaged
  over the window.
- **Where the old term remains:** only in the README's terminology note,
  which says the time-averaged midline level generalizes the classical MESOR
  to a time-varying baseline.
- **Internally,** the key stays "mesor", so no estimator changed.
- **Guard test:** `tests/no_mesor_app_test.R` renders every cosinor output,
  the report and the downloads under both approaches, and fails if the word
  is visible anywhere. A menu's option values are not visible text.

**The (i) buttons.** Every reported quantity has an (i) that opens the
formula computed, written for the fitted model:
- **Where:** the streamlined view, "By question", the DLMO comparisons tab
  and the functional ANOVA tab.
- **What it shows:** the trend term chosen, the harmonics fitted, the log
  scale when on, and whether the group value is a model contrast (mixed) or
  a mean of participants (two-stage).
- **One dictionary,** `server/08n_helpers_equations.R`, restates the code's
  own formulas (for example `dance_trend_factors`) and is not a second
  implementation. The modal comes from `server/95_equations.R`.
- **Test:** `tests/structures_app_test.R` checks, for every structure and
  approach, that every outcome and test key the tables report has its own
  equation (no fallback text), that the streamlined view shows the table's
  own omnibus statistic and level estimate, and that an (i) click opens.

## 28. The overview and the detailed table carry the same tests; the unsampled intervals name their columns

**A mismatch between the overview and the detailed table.**
- **The mixed model:** it reported one rhythm test, all harmonics jointly.
  The overview looked for one test per harmonic, so it showed "—" for the
  24 h and 12 h rows.
- **Now under mixed-effects:** the model also tests each harmonic's
  (cos, sin) pair on its own, with more than one harmonic. These are
  `dance_traj_tests` blocks `harmonic1`, `harmonic2`, ..., each the same Wald
  F of the cell contrasts restricted to that pair. The outcome table lists
  them as `rhythm1`, `rhythm2`, ..., as two-stage already did.
- **The overview's rhythm rows** are built from the table's rhythm tests:
  the joint test (mixed-effects), then each harmonic's own.
- **Test:** `tests/structures_app_test.R` now checks that every rhythm test
  and every primary level and accumulation test of the table appears, with
  its statistic, in the overview.
- **A label fix:** block labels are now looked up through
  `dance_traj_block_label()`. The report's `[[ ]]` lookup stopped on the new
  names.

**Unsampled intervals, made checkable.** A gap is two consecutive sampling
times, among the columns with at least one observation, further apart than
the threshold. The threshold is the user's, or by default twice the median
interval.
- **Columns named:** each gap now carries the two column names either side.
  The summary's assumption line reads "10:00 (+1d)-15:00 (+1d) (no column
  between 'X' and 'Y')", and hovering over "no samples" on the plot says the
  same.
- **Times with the day:** the times show the day, so 10:00 on the second day
  is not read as 10:00 on the first.
- **Switching it off:** a threshold of 0 switches the marking off.

**(i) on the individual tables.**
- **The parameter summary table:** every row has its equation.
- **The individual-results table:** a row of (i) links, one per column,
  including R², % variance and the per-participant zero-amplitude F test.
- **Test:** the harness checks that every link has an equation.

## 29. One trajectory (no between or within factor): estimates, intervals and tests against zero

With no factor there is nothing to compare between groups. The questions
become: what is the trajectory, how precisely is it estimated, and is there
a rhythm or a trend at all.

**Estimates.**
- **Mixed-effects amplitudes and acrophases:** with no design factor the
  mixed model gave none, because `dance_traj_cell_coefs` refused. One
  trajectory is now one cell, "(all)", via emmeans `~ 1`, so it has its
  amplitudes and acrophases with intervals, as every cell does.
- **Agreement:** mixed-effects and two-stage estimates of the same trajectory
  agree, for example a 24 h amplitude of 1.168 vs 1.163 on the test data.

**Tests against zero.** Both are established procedures:
- **Mixed-effects:** a Wald F test (Satterthwaite or Kenward-Roger df, the
  same `dance_traj_L_test` as every block) of the fixed effects against
  zero. It covers all harmonics jointly, each harmonic's (cos, sin) pair,
  and the trend coefficient.
- **Two-stage:** the population-mean cosinor's own tests (Bingham et al.
  1982):
  - a one-sample Hotelling T² on the participants' coefficients,
    F = (n − p)T²/[p(n − 1)] on p and n − p df, all harmonics jointly and
    per harmonic;
  - for the trend, a one-sample t test of the participants' changes.
- **Where they go:** the outcome table carries them, keyed circadian,
  rhythmK and change, so the overview shows them under "Test against zero",
  each with its (i) equation.

**The residual-diagnostics tab** (tab 5) now lists, for one trajectory:
- every estimate with its 95% interval and notes;
- the tests against zero;
- for mixed-effects, the variance components (the random-effect SDs and the
  residual SD);
- for two-stage, the between-participant spread (SD and range) of each
  quantity.

With factors, it points to tab 6.

**A fix in the exported script.** It computed the two-stage outcome table
only when there were groups, so a single trajectory got none. It now
computes it in both cases.

**Test:** `tests/single_traj_app_test.R` runs both approaches and checks:
- the intervals, and the tests against zero recovering a simulated 24 h
  rhythm and trend while finding no 12 h rhythm;
- the tab 5 panel and the overview;
- the equations;
- agreement between the approaches;
- the script reproducing every estimate and test.

## 30. A time origin before the first observation (e.g. wake time)

**The option.** Under Model Specification, "Time origin (t = 0)" is either
the first observation (the default, unchanged) or a clock time before it.
The origin is the last time that clock time occurs at or before the first
observation, so it lies 0 to 24 h before it. There is one origin for
everyone: the harmonics share one axis, and a per-participant origin would
put the participants' acrophases on different clocks.

**What moves and what does not.** Model time 0 is the origin. Both the trend
and the harmonics are anchored there, so the intercept is the midline at the
origin. That is an extrapolation, since nobody was sampled there, and it is
stated as one in the summary, the report, the assumptions and the script.

Every window quantity is still taken over the sampled window [t₁, T], in
hours since the origin, where t₁ is the first observation:
- the time-averaged fitted response;
- the time-averaged midline level;
- the change, f(T) − f(t₁);
- the average rate, the change over (T − t₁);
- the initial rate f′(t₁), the slope at the first observation.

For the saturating trend, A/τ is the slope at the origin.

**The closed forms**, `dance_trend_factors(trend, T, tau, a = t₁)` and
`dance_harmonic_mean_factors(K, P, T, a)`. With a = 0 they reduce to the
earlier forms.

| trend | mean of f over [a, T] | f(T) − f(a) | f′(a) |
|---|---|---|---|
| linear b·t | b(a + T)/2 | b(T − a) | b |
| log b·ln(1+t) | b{[(T+1)ln(1+T) − T] − [(a+1)ln(1+a) − a]}/(T − a) | b[ln(1+T) − ln(1+a)] | b/(1+a) |
| saturating A(1 − e^{−t/τ}) | A[1 − τ(e^{−a/τ} − e^{−T/τ})/(T − a)] | A(e^{−a/τ} − e^{−T/τ}) | (A/τ)e^{−a/τ} |

The harmonic means are:
- cos: [sin(ωT) − sin(ωa)]/[ω(T − a)];
- sin: [cos(ωa) − cos(ωT)]/[ω(T − a)].

Each closed form was checked against numerical integration
(`stats::integrate`) at a = 0 and a = 2.

For the per-cell saturating model (nlme), the delta-method gradient of the
change with respect to τ is A[(a/τ²)e^{−a/τ} − (T/τ²)e^{−T/τ}].

**Why the fit does not change.** With a linear trend, or a saturating trend
with τ held, a shift of origin by d hours is a reparameterisation. The
fitted curve, τ, every window quantity, the amplitudes and the acrophases
are the same, and only two things change:
- linear: the intercept, b₀′ = b₀ − d·b;
- saturating: A′ = A·e^{d/τ} and b₀′ = b₀ + A − A′.

For mixed-effects models this holds up to lme4's convergence. With a
random-effect variance near zero, the optimiser's stopping point and the
Satterthwaite df depend on the parameterisation. In the test the estimates
agree to about 1e-4 and the F statistics to 1.4%.

So for a rising trend, the extrapolated level at wake is **lower** than at
the first observation, and A is larger. It is higher only for a falling
trend. τ itself does not change with the origin, but the share of A reached
by a given hour does. A logarithmic trend, ln(1 + t), changes shape with its
origin, so its fit changes too; the panel says so.

**Script.** The exported script:
- states the origin;
- carries time_vec on the origin's axis;
- passes `t_origin = 0` to `dance_ts_fit_all()` and `dance_ts_model_set()`
  (the model set now also gets its `k_max`);
- uses `t0 = 0` for the mixed model, as before.

Also in this commit, the exported script's per-participant table
(`cosinor_params`) and its bootstrap print call the fit's constant
`intercept`. They had called it `midline_level` since the MESOR rename.

**Test:** `tests/time_origin_app_test.R` runs mixed-effects and two-stage,
each with a linear and a saturating trend, with an origin at 06:00 against a
recording starting at 08:00. It checks:
- that every window quantity and its group test is unchanged;
- that the intercept follows the reparameterisation above;
- that the assumptions and the equations state the origin;
- that the script reproduces every estimate;
- that an invalid time is refused.

**Drawing the curves from the time origin** (display only; no estimate,
test or script number changes). With an origin set, one setting beside
it, "Draw the fitted curves from the origin", draws every fitted-curve plot
from the origin. That covers the Estimated trajectories views (full
trajectory, trend and rhythm) and the Fitted curves tab. It is ticked by
default and can be switched off.

When drawn from the origin:
- **Trend view:** f(t) − f(0), zero at the origin, not at the first
  observation.
- **Extrapolated interval:** the interval from the origin to the first
  observation is shaded, labelled "before the first sample", and the curves
  are dotted over it.
- **Annotations:** they stay on the sampled window and start on the drawn
  curve at the first observation.
  - The total-change line runs from f(t₁) to f(T).
  - The tangent's slope is the rate at t₁.
  - The asymptote is drawn at A − f(ref), and the 63.2% point at t = τ
    after the origin.

**A fix:** the mixed-effects plot had passed the origin, not the first
observation, as the annotations' window start. That put the average-rate
line's start at the origin.

**Test:** in `tests/time_origin_app_test.R`, for both approaches:
- the trend is 0 at the origin;
- the change line starts on the curve at the first observation;
- the interval before it is marked;
- unticked, curve and line start at 0 at the first observation.

## 31. Independent against continuous-time AR(1) residuals: an audit of the pair

**What was required.** The two residual options of the mixed-effects cosinor
must use the same observations, the same fixed-effect design matrix and the
same random structure, so that they differ only in the residual covariance.
Every comparison of the two (likelihood ratio, AIC, BIC) must come from one
pair of fitted models.

**How the AR(1) model is fitted (unchanged).** For a given φ, the CAR(1)
model is the ordinary mixed model fitted to the whitened data W·y, W·X and
W·Z (Prais–Winsten for unequal gaps; W R W′ = I). Its log-likelihood is the
whitened one plus log|det W|.
- **ML:** the identity is immediate.
- **REML:** the restricted term log|X′V⁻¹X| is invariant under the
  transformation (X* = WX, V* = WVW′), and log|V*| = log|V| + 2 log|det W|.
  So L_R(y) = L_R*(Wy) + log|det W|, in lme4's convention (no log|X′X| term)
  for both models alike.
- **φ:** profiled on a 20-point grid in −log(1 − φ), including φ = 0 exactly,
  then refined between the neighbours of the best grid point.

**The audit** (`tests/testthat/test-traj-car1-pair.R`; simulated data with
CAR(1) residuals, uneven times and missing values, with and without a
within-participant factor, five trend set-ups). What already held:
- **Same model parts:** the same rows, the same unwhitened X (and the AR(1)
  fit's X is exactly W·X), the same random structure.
- **φ = 0:** at identical θ, the whitened REML deviance at φ = 0 (and at
  φ = 1e−8) equals the independent model's to 1e−8.
- **REML correctness:** both REML log-likelihoods equal a dense direct
  computation, −½[(n − p) log 2π + log|V| + log|X′V⁻¹X| + r′V⁻¹r], to 1e−7.
- **Optimum:** the AR(1) optimum is above the independent one. Started from
  the independent estimates at φ = 0.001, and at φ̂, the optimizer reaches the
  same optimum as from the default start, to 1e−6.
- **Parameter count:** φ is counted as one parameter more.
- **τ held or no trend:** the reported likelihood ratio equals 2·Δℓ of the two
  reported fits.

**The inconsistencies found and corrected:**
1. **The random structure could change between the two options.** The final
   REML refit at φ̂ was not capped at the independent model's
   random-effects rung (`min_rung` without `max_rung`). If that structure
   failed there, the fit stepped down the ladder. The AR(1) model then had a
   different random structure from the independent one, and its likelihood
   ratio, computed on the profile's structure, belonged to neither reported
   fit. Shown by fault injection: rung 2 was reported against rung 1, with a
   likelihood ratio of 54.68 against 51.04 from the reported pair.

   The joint τ–φ rounds re-chose a structure in each τ profile in the same
   way. Both now use the independent model's structure only. If it fails, the
   AR(1) fit is refused with a message; a different structure is never
   reported.
2. **A φ estimate bounded by non-convergence was reported as an estimate.**
   Grid points whose fit did not converge were dropped. The best remaining
   point could then be the edge of the converged region, and it was reported
   as "at the top of the search (0.98)" even when that edge was 0.186. A
   best point next to a failed grid value is now flagged
   (`bounded_by_failure`) and reported as a bound set by the optimizer.
   "At the top" now means the grid's last value.
3. **With τ estimated, the likelihood ratio and the AIC/BIC came from
   different models.** The likelihood ratio was ML (each model at its own best
   τ), while the printed AIC/BIC were the REML values of fits whose τ, and so
   whose fixed-effect design, can differ. Those REML values are not
   comparable. There is now one comparison table (`residual_ar$comparison`),
   and the likelihood ratio, AIC and BIC all come from its two rows:
   - with τ held, or no saturating trend: the two reported fits, by REML;
   - with τ estimated: one ML refit at each model's estimate, on the same
     random structure.

   The table is printed in the summary, the report and the exported script.
4. **An estimated τ was not counted as a parameter** in the fitted model's
   AIC/BIC, though the model comparison counted it. It is counted now; a held
   τ is not. A fit whose maximized likelihood is below the independent
   model's, which it contains, would be flagged as not having reached the
   maximum. None of the audited cases is.

**What is unchanged:** the model specification, the estimators and every
estimate where none of the four cases applies. In the audited cases, the
reported fits, φ̂ and the REML likelihood ratios are identical before and after.

## 32. Per-participant saturating fits failed whenever the response fell

**The report.** Pilot data, DOZEPROB, a saturating trend with 24 h and 12 h
harmonics, mixed-effects: "9 of 107 subjects failed to fit (need 8+ valid
points)". All nine had all 13 points.

**The cause.** In the per-participant nonlinear fit, the start value for A
is `y_range * 1.5` when the participant's linear slope is positive. Otherwise
it is `lin_slope * t_max`, and `lin_slope` was `coef(lm(...))[2]`, which
keeps the regressor's name. nls names each parameter after its start value,
so the parameter became `A_sat.t_shifted`, and every reader looked up
`A_sat`.
- **τ held:** the fit was discarded ("Result extraction failed: missing value
  where TRUE/FALSE needed", from `append(after = match("A_sat", ...))`).
- **τ estimated:** the fit was kept with A = NA, so that participant's trend
  change, average rate and time-averaged midline level were lost.

So every participant whose response fell over the recording (a non-positive
linear slope) lost the saturating trend. The nine in the pilot are exactly
its nine participants with a non-positive slope.

**Who was affected:**
- **Two-stage approach:** those participants were excluded from the group
  results (τ held), or dropped out of the trend summaries (τ estimated).
- **Mixed-effects approach:** the model is fitted to every participant's data,
  so it was unaffected. Only the per-participant views lost them.
- **The two-stage model comparison:** its candidate fits use the same
  function.
- **The exported script:** it calls the same function. The fix reaches it.

**The fix.** The coefficient is unnamed, and every start value is stripped of
names before fitting. The warning now gives the reason for each failure
(`dance_ts_failure_text()`):
- too few observed points;
- no variation (every value equal);
- the fit did not converge.

Under mixed-effects it also says that only the per-participant views are
missing those participants. In the pilot, the one remaining failure, S70, has
no variation at all (15 at every time point).

**Test:** `tests/declining_trend_app_test.R`. Half the simulated participants
fall, and one is flat. Under τ held and τ estimated, it checks that:
- every participant with variation is fitted, and the falling ones get a
  negative A;
- the flat one is reported as having no variation;
- the exported script fits the same participants and reproduces every
  estimate.

On the code before the fix, the test fails: 11 of the 12 falling participants
are not fitted with τ held.

## 33. Change since the time origin; plot and interface changes

**Change since the time origin** (`change_origin`, an estimate). When an
origin before the first observation is set, the detailed statistics also
give the trend's change from the origin to the last observation, f(T) − f(0).
It is flagged as extrapolated before the first observation. The window
change, f(T) − f(t₁), stays the primary quantity.
- **Mixed, τ shared or held (or a linear trend):** it is a linear functional
  of the same trend coefficient, the window change times
  f(T)/(f(T) − f(t₁)). Its Wald test is the same test (same F, same p).
- **Mixed, τ per cell (nlme):** A(1 − e^(−T/τ)), with the delta-method
  gradient (1 − e^(−T/τ), −A(T/τ²)e^(−T/τ)·dτ/dx) and its own Wald test.
- **Two-stage:** each participant's own trend(T) − trend(origin), then
  group means and the same one-way ANOVA. With τ held, or a linear trend, the
  values are proportional to the window change, so the F is identical.
- **Tests:** `tests/time_origin_app_test.R` checks that the row is absent
  without an origin, and that with one it equals the window change × the
  factor (1.0667 linear, 1.1740 saturating), with the same test, for both
  approaches. The script round-trip covers the new rows.

**Plots and interface** (display only; no estimate changes):
- **Fitted curves, one participant:** their raw data are their own legend
  entry, so hiding the group's curve no longer hides them.
- **Raw points:** larger (size 8, outlined).
- **Components:** each harmonic and the trend can now be shown under the
  mixed-effects approach too, from the model's own component predictions,
  each its own legend entry.
- **Bands:** the Estimated trajectories plot has a band option "None".
- **Group Structure:** the plot now shows the raw curves with each group's
  mean on the clock axis. The old plot put the curves on a 0–1 axis and the
  means on clock hours, so the means fell off the edge. The duplicate plot in
  Data Preview is removed.
- **Participant table:** the note on the mixed model's per-participant values
  is rewritten in plain language.

**Explanations behind an (i).** Every explanatory note is now hidden behind
an (i) that opens a pop-up (`dance_help()`, `server/00a_help.R`). That means
every `helpText()` in the app, and the long grey notes of the cosinor tab.
The browser moves each (i) next to the label of the control just above it
when both are in the same container. Messages that must be seen (nothing
selected yet, a missing package, a test's method, warnings about the time
axis) stay visible through `dance_status()`.
- **Tests:** `tests/testthat/test-help-popups.R` and
  `tests/fitted_curves_view_test.R`.

**Plot layering (display only).** Raw data are drawn large and outlined for
one participant, and small and faint (size 5, 22% opacity) when many curves
are shown. Every cosinor curve plot is drawn back to front, in this order:
1. the data (raw points, every participant's own curve);
2. the bands;
3. the fitted lines (`dance_plot_layers()`).

So the fits are never buried under the data.

## 34. One sign for every pairwise row; publication-ready text in the report

**The sign of a pairwise difference.** A row labelled "a vs b" is now
a − b for every quantity, under both approaches.
- **The bug.** Under the mixed model, the amplitude and acrophase rows were
  b − a, while the level and change rows were a − b. With Group 2's
  amplitude larger than Group 1's, the row "1 vs 2" showed a positive
  difference. The estimates were right; only the sign of those rows was
  reversed.
- **Where it came from.** The low-level helpers `dance_traj_phase_contrast()`
  and `dance_traj_pair_joint()` (`server/08f`) return second − first. Their
  own tests, including the GLMMcosinor crosscheck, pin that.
- **The fix.** `dance_traj_contrasts()` (`server/08g`) now flips the
  estimate, the statistic and the CI of those rows, so the table it hands on
  is a − b. The p values do not change. The two-stage pairwise rows
  (`dance_ts_pairwise()`) were already a − b.
- **Two-stage amplitude pairwise.** Under the two-stage approach, the
  amplitude pairwise row compares the mean of the participants' own
  amplitudes (Welch). The group amplitude in the outcome table is the
  population-mean cosinor's amplitude (the vector mean, Bingham et al.,
  1982). The two estimators differ, so the pairwise difference need not
  equal the difference of the two table values. The same holds for
  acrophase (Watson-Williams on the participants' acrophases).
- **Tests:**
  - `tests/pairwise_sign_app_test.R` checks that "a vs b" = est(a) − est(b)
    for level, change, amplitude and acrophase under both approaches;
  - `tests/glmmcosinor_crosscheck_test.R` checks the amplitude contrast
    against GLMMcosinor with the sign that "a vs b" now carries.

**Publication-ready text (report only, no new number).**
`dance_pub_cosinor()` (`server/08p_helpers_pubtext.R`) writes a "Statistical
Analysis" section, a "Results" section and a "Table 1" with a Note into the
analysis report, after the results tables.
- **Its numbers.** Every number is read from the outcome table
  (`harmonic_outcomes()`) and the pairwise table
  (`harmonic_pairwise_all()`). Nothing is computed anew, so the exported
  script needs no change: it already reproduces both tables.
- **Its wording, by rule:**
  - a test is called significant when p < alpha (.05);
  - a significant pairwise difference is written positive, with the larger
    group named first; the others are listed with their p;
  - amplitude pairwise comparisons are called exploratory, with a caution
    when the joint rhythm test is not significant or nothing is adjusted;
  - the summary names the domains whose tests were significant.
- **Its caveats.** The report says the text is for adaptation and must be
  read against the results. It also says the components describe the
  observed trajectories, not endogenous circadian or homeostatic processes.
- **Tests:** `tests/pub_text_app_test.R` checks:
  - the section exists under both approaches and for numeric group levels;
  - the level estimates and CIs equal the outcome table's;
  - the level contrast is stated positive, with the larger group first;
  - the rhythm joint test is stated with its F;
  - Table 1 uses the group labels;
  - every sentence starts with a capital, with no doubled nouns.

## 35. Participants without a group; a sparse R-squared; labels and panel grid

**Participants with no group value (result-relevant).**
- **What happened.** A participant with a blank or missing group was kept
  in the mixed model's data frame. lme4 left those rows out of the fit on
  its own (`na.omit`), but everything built from the frame still carried
  them. So:
  - the residuals were longer than the fitted values ("longer object length
    is not a multiple");
  - a blank cell was read as a group called "", which has no colour, and
    the Fitted curves plot stopped ("subscript out of bounds").
- **Now.**
  - `dance_traj_long()` reads a blank label as missing.
  - `dance_traj_spec()` drops participants with a missing value for a
    factor of the design, once, and records it in `spec$missing_note`
    (how many, which factor, which participants).
  - The note is printed in "What this run assumed", so it reaches the
    summary, the report and the exported script.
  - Two-stage (`dance_ts_frame()`) also reads a blank label as missing, as
    the population cosinor already did (`dance_group_audit()`).
- **Effect on results.** The mixed model's estimates do not change: lme4
  had already left those rows out. Two-stage group summaries no longer
  include a group "".
- **Tests:** `tests/testthat/test-traj-missing-factor.R`.

**R-squared with a sparse Z (same formula).** `dance_traj_r2()` made the
random-effects design matrix Z dense. Z is observations × (participants ×
terms), so on a large data set this needed 1.3 GiB. The sum of z_i' Σ z_i is
now computed with Z kept sparse. A row has non-zeros under one participant
only, so summing the term-by-term products over all participants gives that
participant's product. The test checks the result against the dense sum to
1e-10, and `test-traj-r2.R` still passes.

**Display only (no script change).**
- The pairwise "Compare" menu showed "NA" for the change since the time
  origin, which had no label. It now has one, and a name without a label is
  shown as itself.
- The panel titles under "Intercept, midline level and change" now name
  `fitted_mean`, `change_origin`, `avg_rate` and `rate0` in words.
- Dot-and-whisker panels wrap after 4 per row, and the plot grows by one
  row height per row.
- With more than 200 participants, the "Show" menu of the Fitted curves
  loads its options from the server (server-side selectize).
