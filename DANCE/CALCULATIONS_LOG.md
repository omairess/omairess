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
