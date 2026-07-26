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
