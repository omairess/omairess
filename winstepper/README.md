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
  with Rasch-Thurstone (2.3) and modal (2.1) variants.
- **DGF** (Advanced ▸ DGF): Differential Group Functioning — the interaction of
  the item scales with a person classification, alongside item-level DIF.

## Files

- `app.R` — the WINSTEPPER Shiny UI/server (bslib).
- `rasch_engine.R` — the Rasch engine (no Shiny dependency; reused unchanged).
- `winsteps_plots.R` — all figures (base graphics; reused unchanged).
- `house_modules.R` — shared house-style modules (reused unchanged).
- `CLAUDE.md` — project memory and non-negotiable conventions.

## Honest scope

Not byte-identical to WINSTEPS 5.11 (agreement ≈ 2 decimals on well-behaved
data). No `STBIAS=` bias correction; no anchoring, `CUTLO/CUTHI`, subset
detection, keyform/scalogram tables, DPF, or non-uniform DIF. Independent
software; not affiliated with, endorsed by, or derived from WINSTEPS®, which is
John M. Linacre's software. Cite the engine you actually ran.
