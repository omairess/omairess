# CoDA sleep diaries — a Shiny front end for `multilevelcoda`

Bayesian multilevel compositional data analysis of sleep diaries with
[multilevelcoda](https://florale.github.io/multilevelcoda/) (Le, Stanford,
Dumuid & Wiley, 2025). Built to the shiny-house-style rules: every step is
recorded, and the app exports a stand-alone R script that reproduces the
on-screen numbers exactly.

## Run it

```r
# once: a C++ toolchain for Stan (macOS: run `xcode-select --install` in Terminal)
install.packages("multilevelcoda")      # pulls in brms, compositions, ...
shiny::runApp("multilevelcoda-app")     # from the repository root
```

Missing packages are installed on first start; installed ones are never
updated by the app. Tested with R 4.3, multilevelcoda 1.3.3, brms 2.20,
rstan 2.32, shiny 1.14, shinyWidgets 0.9.

## Workflow (one tab each)

1. **Data** — upload one or several files. Each sheet or table is read in one
   of three layouts, and all files are stacked into one row per participant-night.
2. **Diary & composition** — map columns to diary items (guessed from Dutch,
   French and English names), derive the sleep parameters, choose the
   composition, handle zeros, pick outcomes and align them in time.
3. **Model** — `complr()` + `brmcoda()`: the composition predicts an outcome
   (between- and within-person ilrs) or is itself the outcome.
4. **Substitution** — `substitution()`: expected change in the outcome when
   minutes move between parts.
5. **Pipeline & export** — pipeline log, reproducible script, data (.csv /
   .xlsx) and all R objects (.rds).

## Supported inputs

| Your file | Layout to choose | ID from |
|---|---|---|
| Study workbook, one row per night, a participant column (e.g. *SVT slaapdagboeken*) | One row per night | A column |
| One person per row, day-numbered columns (`QS_1 … QS_14`) | One row per person | A column |
| CBT-I diary per patient: items in rows (*Heure de coucher*, *QS*, *QR*, *FJ*), days in columns; weeks as separate sheets or stacked in one sheet (e.g. *MARISKA.xlsx*) | Items in rows, days in columns + **All sheets** | File name — upload all patients of a group at once to pool them |
| Sleep-diary app export, `.json` (e.g. *Justine.json*) | One row per night | File name |

Also `.csv/.txt` (delimiter detected), SPSS `.sav`, `.rds`, `.RData`. Title
rows above the header are detected; accents are removed from column names so
they can enter model formulas (`Qualité du sommeil` → `Qualite_du_sommeil`).

For JSON, every array of records (or object of records keyed by date) is
offered; the largest is preselected. Top-level fields such as a participant
code are kept as `meta_*` columns.

Times are read from `23:30`, `23h30`, `23u30`, `11:30 pm`, `2330`, Excel
time cells and Excel date-times. Durations are read from minutes,
`1h30`, `0:45`, `20 min` and ranges such as `10-15` (midpoint). Whether
`23.30` means 23:30 or 23 h 18 min is a setting.

## Sleep parameters (Consensus Sleep Diary; Carney et al., 2012)

- TIB = out of bed − into bed (across midnight)
- TST = (final awakening − lights off) − SOL − WASO, or the recorded TST
- EMA = out of bed − final awakening; PRE = lights off − into bed
- SE = TST / TIB × 100
- A night is excluded (and listed with the reason) when TIB or TST cannot
  be computed, TST ≤ 0, TST > TIB, TIB ≥ 24 h, or the clocks are out of order.

## Compositions

| Option | Parts | Closed to |
|---|---|---|
| 24 h, 3 parts | Sleep, WakeInBed (TIB − TST), OutOfBed (1440 − TIB) | 1440 min |
| 24 h, 5 parts | TST, SOL, WASO, EMA (+ pre-lights-off), OutOfBed | 1440 min |
| Time in bed, 4 parts | TST, SOL, WASO, EMA | mean TIB |
| Custom | any numeric columns (e.g. accelerometry minutes) | your total |

Log-ratios need strictly positive parts. Zeros (frequent for WASO and SOL)
are replaced multiplicatively with 0.65 × a detection limit you set
(Martín-Fernández et al., 2003), or the night is dropped; the count is
logged. Many replaced zeros make a part tiny at the mean composition, which
limits how many minutes can be moved out of it in substitution; the
Substitution tab shows the compositional means for that reason.

## Methodological notes

- **Between vs within.** With several participants the composition enters
  as between-person ilrs (a person's usual composition) and within-person
  ilrs (a night's deviation from it), each with its own coefficients.
  With fewer than about 10 participants the between-person part is very
  uncertain; the app says so.
- **One patient.** A single diary gives a single-level model of the nights
  (no random effects; aggregate-level substitution). Nights are serially
  dependent, so these estimates describe that patient. Pool a CBT-I group
  by uploading all patients' files together.
- **Timing of ratings.** QS and QR in a morning diary describe the night on
  the same row. If FJ (form of the day) refers to the day after that night but
  sits on the next row or column, select it under *Take from the NEXT row*.
- **Priors and convergence.** brms default priors are used and printed in
  the script (`prior_summary`). The Model tab reports R-hat, ESS ratio and
  divergent transitions. Raise adapt_delta or the iterations if they are poor.
- **Interpretation.** Substitution estimates are model-based contrasts
  around the reference composition, not causal effects of changing
  behaviour.

## Files

```
multilevelcoda-app/
  app.R                         navbar UI + server wiring
  R/00_packages.R               package guard + version ledger
  R/01_recorder.R               step recorder -> pipeline log + R script
  R/02_diary_helpers.R          import, diary derivation, composition, figures
                                (plain R; embedded verbatim in the exported script)
  R/03_data_module.R            Data tab
  R/04_plot_export.R            colour pickers, PNG/PDF/SVG export
  R/05_diary_module.R           Diary & composition tab
  R/06_model_module.R           Model tab (complr + brmcoda)
  R/07_substitution_module.R    Substitution tab
  R/08_exports_module.R         Pipeline & export tab
```

## References

- Le, F., Stanford, T. E., Dumuid, D., & Wiley, J. F. (2025). Bayesian multilevel compositional data analysis: Introduction, evaluation, and application. *Psychological Methods*. https://doi.org/10.1037/met0000750
- Le, F., Dumuid, D., Stanford, T. E., & Wiley, J. F. (2025). Bayesian multilevel compositional data analysis with the R package multilevelcoda. *Multivariate Behavioral Research*. https://doi.org/10.1080/00273171.2025.2565598
- Carney, C. E., et al. (2012). The Consensus Sleep Diary: Standardizing prospective sleep self-monitoring. *Sleep, 35*(2), 287–302. https://doi.org/10.5665/sleep.1642
- Aitchison, J. (1982). The statistical analysis of compositional data. *Journal of the Royal Statistical Society: Series B, 44*(2), 139–177.
- Egozcue, J. J., Pawlowsky-Glahn, V., Mateu-Figueras, G., & Barceló-Vidal, C. (2003). Isometric logratio transformations for compositional data analysis. *Mathematical Geology, 35*(3), 279–300.
- Martín-Fernández, J. A., Barceló-Vidal, C., & Pawlowsky-Glahn, V. (2003). Dealing with zeros and missing values in compositional data sets using nonparametric imputation. *Mathematical Geology, 35*(3), 253–278.
- Dumuid, D., et al. (2019). The compositional isotemporal substitution model. *Statistical Methods in Medical Research, 28*(3), 846–857.
- Bürkner, P.-C. (2017). brms: An R package for Bayesian multilevel models using Stan. *Journal of Statistical Software, 80*(1), 1–28. https://doi.org/10.18637/jss.v080.i01
