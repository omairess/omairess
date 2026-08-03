# WINSTEPPER

A Shiny app for WINSTEPS-style Rasch measurement — a redesigned front end built
on the audited R-Winsteps estimation engine. Dichotomous Rasch, Andrich Rating
Scale and Masters Partial Credit models via JMLE (UCON), with fit statistics,
separation/reliability, rating-scale category structure, Wright maps, a PCA of
residuals (Table 23), DIF (Table 30), and a fully reproducible R-script export.

This folder is self-contained and independent of the other apps in the repo.

## Run

```r
# from the repository root
shiny::runApp("winstepper")
```

The first launch installs any missing packages (it never auto-updates existing
ones). Then either upload a wide data file (`.csv/.txt/.xls/.xlsx/.sav/.rds/.RData`)
or click **Load built-in demo data** on the Data tab, select your item columns
and a model on the Estimate tab, and press **Estimate**.

### Valid category codes (`CODES=`) and recoding (`NEWSCORE=`)

On the Estimate tab you can list the response codes that count as valid
categories, exactly like WINSTEPS `CODES=` (e.g. `0,1,2,3`). Any value **not**
in the list is treated as missing / not administered before estimation. Leave
it blank to keep every observed value. The sidebar shows a live count of how
many responses fall outside the codes.

`NEWSCORE=` positionally recodes the valid codes — e.g. `CODES=0,1,2,3` with
`NEWSCORE=0,1,1,2` collapses categories 1 and 2. Both are written into the
exported reproducible R script.

### Grouped / mixed model — item scales

Choosing the **Grouped / mixed** model reveals a per-scale editor: set the number
of scales, then for each scale pick which items belong to it and (optionally) give
that scale its own `CODES` / `NEWSCORE`. This is how you assign items to shared
rating-scale structures (WINSTEPS `ISGROUPS=`) while others stay partial-credit.

### Keyform (Table 2.2) and DGF (Table 33)

- **Keyform** (Results ▸ Keyform): the general expected-score keyform (Table 2.2),
  with Rasch-Thurstone (2.3) and modal (2.1) variants, and an optional person
  histogram underneath sharing the same logit axis.
- **DGF** (Advanced ▸ DGF): Differential Group Functioning — the interaction of
  the item scales with a person classification, alongside item-level DIF.

### Dropping misfitting persons

On the Persons tab you can flag persons by a fit statistic (Outfit/Infit
MNSQ or ZSTD above a cut), edit the exclusion list by hand, then press
**Re-estimate without excluded persons**. The model is recalibrated from
scratch without them, so item difficulties, thresholds and every downstream
table update. The exclusions are written into the exported script.

Report this honestly: removing misfitting persons improves fit almost by
construction. State how many were removed, on what criterion, and what changed.

### Person–item barchart

On the **Wright map** tab, beneath the Wright map itself: person-measure
histograms — the whole sample plus one panel per level of a classification you
choose (sex, group, …) — beside every item drawn as a vertical string of
numbered thresholds, all sharing one logit axis, so a person's location reads
straight across into the item thresholds. It mirrors the WINSTEPS Plots options:
measure-relative values (Andrich, Thurstonian, half-point, or
maximum-probability full-point thresholds), sort order (entry / measure / alpha,
ascending or descending), right-side size and title, and whether to mark each
item's own measure.

### Suggested rescore

When a rating scale misbehaves, the **Rating scale** tab shows a *Suggested rescore*
card: a proposed collapse in the same comma-separated format the rescore box expects,
for example

```
CODES:    0,1,2,3,4,5,6,7,8,9,10
NEWSCORE: 0,1,1,1,1,2,2,2,3,3,3
```

together with the resulting bands, their sizes, and why it proposes them — plus an
alternative one category coarser.

**What it optimises.** It keeps **as many categories as possible** whose thresholds are
still at least `min separation` logits apart (default 1.0), preferring evenly spaced
bands. That is deliberately the opposite of "merge until the threshold-advance guideline
is satisfied": on a scale whose categories are genuinely compressed, no collapse can ever
satisfy that guideline, so such a rule collapses everything into two poles and one huge
middle — which destroys the information function and person separation. Raise the
separation for a coarser scale, lower it to retain more categories.

The basis is the observed cumulative thresholds, `-logit(P(X ≥ k))`. Under an adjacent
merge these are exactly a subset of the originals, so every candidate collapse is scored
without any re-estimation, and the search over all 2^m partitions is exhaustive rather
than heuristic.

**It is only a suggestion.** Nothing is applied and nothing is re-estimated. The button
copies the strings into the CODES / NEWSCORE boxes on the Estimate tab and takes you
there; edit them however you like, then press **Estimate** yourself. Since thresholds are
re-estimated once categories change, the intended loop is suggest → apply → re-estimate →
re-read this tab.

There is no AI service behind this: it is plain R, works offline, and costs nothing.

### Threshold advance criterion

The category-quality check uses the **category-count-dependent** minimum
threshold advance, `ln((k+1)(m+1−k) / (k(m−k)))`, not a flat 1.4 logits — 1.4 is
only the three-category case (2 ln 2 = 1.386). For 4 categories the requirement
is 1.10/1.10, for 5 it is 0.98/0.81/0.98, and so on (Linacre, *Rasch Measurement
Transactions* 2006, 20:1, p. 1052).

### CMLE alongside JMLE

Estimate ▸ **Estimation method** offers conditional maximum likelihood as well
as the default JMLE. CML conditions each response pattern on its person raw
score, which is sufficient for θ under the Rasch model, so the person parameters
drop out of the likelihood entirely. The item estimates are then consistent as
N grows with the test length fixed — JMLE's ~L/(L−1) spread inflation, which
this project has no `STBIAS=` correction for, is absent by construction.
Persons are measured afterwards by ML with the items anchored, so every fit
statistic, reliability and residual table downstream is unchanged.

Covers dichotomous, RSM and PCM data, and mixed item groups. Slower than JMLE,
and proportionally slower again when the data contain many distinct
missing-data patterns — each pattern conditions on its own item set and so
needs its own elementary symmetric functions.

Verified by `test_cmle.R`: the symmetric functions match brute-force
enumeration of every response pattern to 1e-10; a two-item test matches the
closed form; the analytic gradient matches a numerically differentiated
brute-force likelihood; the CML equations (observed = conditionally expected
sufficient statistic) hold at the solution; parameters are recovered from
simulated data; and the JMLE/CML spread ratio comes out near L/(L−1).

### Labels always fit

Figures size their margins from the labels they actually have to draw, at the
label size currently set on the Figure style tab, rather than from a fixed
number of margin lines. Text that cannot fit even a capped share of the figure
is truncated with an ellipsis rather than clipped. `test_collapse.R` section 8
checks the width estimate against `strwidth()` and fails if it ever
under-shoots.

### Interface font size

Settings ▸ Figure style ▸ **App appearance** scales the whole interface
(root `html` font size; Bootstrap sizing is rem-based) and the table/console
font separately. Figure text has its own size sliders.

## Files

- `app.R` — the WINSTEPPER Shiny UI/server (bslib).
- `rasch_engine.R` — the JMLE engine (no Shiny dependency; reused unchanged).
- `winstepper_cmle.R` — conditional maximum likelihood (no Shiny dependency).
- `winsteps_plots.R` — all figures (base graphics).
- `winstepper_extras.R` — keyform, DGF, person-item barchart, engine overrides.
- `house_modules.R` — shared house-style modules (reused unchanged).
- `test_cmle.R`, `test_collapse.R` — the test suites; run both before shipping.
- `CLAUDE.md` — project memory and non-negotiable conventions.

## Honest scope

Not byte-identical to WINSTEPS 5.11 (agreement ≈ 2 decimals on well-behaved
data). No `STBIAS=` bias correction for JMLE — use CMLE if that bias matters,
though our CML is an independent implementation and will not match WINSTEPS'
to the last decimal either. No anchoring, `CUTLO/CUTHI`, subset
detection, keyform/scalogram tables, DPF, or non-uniform DIF. Independent
software; not affiliated with, endorsed by, or derived from WINSTEPS®, which is
John M. Linacre's software. Cite the engine you actually ran.
