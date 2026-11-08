# BioMas 1.3.5 — Biomathematical Sleep Modeling

Version 1.3.5 adds **time-in-bed (bedtime-interval) schedule input** from actigraphy or sleep diaries and makes the model's **initial conditions explicit as a baseline time in bed**, defaulting to the baseline used in the 2024 model paper. The model equations, parameters and calculation functions are unchanged from 1.3.4.

## Why time in bed, not total sleep time

The 2024 model (McCauley et al., *Front. Environ. Health* 2024;3:1362755, [doi:10.3389/fenvh.2024.1362755](https://doi.org/10.3389/fenvh.2024.1362755)) was calibrated on laboratory **sleep opportunities**: scheduled lights-out to lights-on intervals, such as "10 h sleep opportunities (22:00–08:00)". The model has only two states, asleep and awake, and treats the whole opportunity as sleep. Time awake in bed is absorbed into the calibrated parameters. Entering total sleep time (TST), or splitting nights at scored wake epochs, therefore departs from the calibration convention and usually predicts more impairment.

## Initial conditions: baseline time in bed

*Parameters → Initial conditions: baseline time in bed*

- The fields are now **Habitual bedtime, lights out** and **Habitual rise time, lights on**. Enter time in bed, not TST.
- The default is now **22:00–08:00 (10 h)**, the baseline of studies Cal1, Cal3 and Val1–Val3 in the 2024 paper (Cal2 used 23:30–07:30). The previous default was 23:00–07:00 (8 h).
- **Reset to paper default (22:00–08:00)** restores it. A line below the fields shows hours in bed and awake per baseline day.
- The app repeats this baseline day until it reaches a stable state, then uses that state as p, u, h and κ at study start. This remains an assumption about pre-study sleep, not a measured individual state.

Changing the baseline changes the predictions. In the paper's Cal1 protocol, the 2024 model predicts a mean KSS of 3.31 between 10:00 and 20:00 on the first study day with a 22:00–08:00 baseline, and 3.78 with 23:00–07:00. Choose the baseline from your participants' actual pre-study time in bed when it is known.

## Bedtime-interval input (actigraphy rest intervals or diaries)

*Schedule → Input mode → Actigraphy / diary bedtime intervals (time in bed)*

Upload a CSV or Excel file with **one row per rest interval**: bedtime (lights out) and rise time (lights on). Typical sources are actigraphy rest-interval exports (for example Actiware `Start Date`/`Start Time`/`End Date`/`End Time` with `Interval Type`, or ActiLife `In Bed Time`/`Out Bed Time`) and sleep diaries. An example file is included: `examples/bedtime_intervals_example.csv`.

1. **Time format**
   - *Local date-time* (most exports): every row gets the UTC offset of the Study start. Turn on **Use dated schedule / actigraphy** and enter the Study start with its offset. Choose the **date order** for non-ISO dates such as `18/10/2026`; day/month order is never guessed. 12-hour clocks with AM/PM are accepted.
   - *ISO date-time with UTC offset*: each row carries its own offset.
   - *Elapsed hours from t = 0*.
2. **Columns** are proposed from the header names and can be changed. Dates and times may be in one column or in separate date and time columns. Excel date and time cells are read without conversion.
3. **Interval type**: if the file has an interval-type column, keep only **REST** rows (pre-selected when present). Sleep-scored intervals would give TST rather than time in bed. Mixing both types produces overlapping intervals, which are rejected.
4. **Per-night audit table**: every interval is listed with its elapsed and clock times, time in bed, and role:
   - *Pre-study (baseline)*: ends at or before t = 0. Not simulated.
   - *Study schedule*: simulated whole as a sleep opportunity.
   - *Study schedule (clipped)*: spans t = 0 or the study end and is cut at that boundary.
   - *After study end (ignored)*.
5. **Use pre-study nights as baseline time in bed** sets the initial-conditions fields to the circular-mean bedtime and rise time of the pre-study nights, rounded to 5 minutes. Circular means handle nights on either side of midnight. The status line also shows the circular SD of bedtimes and rise times, so irregular sleepers are visible.

Validation is strict: a row with a missing bedtime or rise time, a rise time before its bedtime, time in bed over 24 h, or overlapping intervals stops the run and names the file rows. Blank rows are skipped. Runs record the column mapping, filter, format, per-night table and baseline assumption in the exported bundle.

The epoch-by-epoch mode is renamed **Actigraphy epochs (sleep/wake scoring, experimental)**. It is unchanged, and its help text now points to the bedtime-interval mode for calibrated input.

## Run

Stop the previous app. Extract the ZIP and run, from its parent folder:

```r
shiny::runApp("BioMas_1.3.5")
```

or open `BioMas_1.3.5/app.R` in RStudio and click **Run App**. The dashboard header should show **BioMas 1.3.5**.

Required packages:

```r
install.packages(c(
  "shiny", "shinydashboard", "shinyWidgets", "DT", "plotly", "ggplot2",
  "dplyr", "tidyr", "readr", "readxl", "haven", "lubridate",
  "RColorBrewer", "deSolve", "digest"
))
```

## Tests

From the extracted `BioMas_1.3.5` folder:

```r
install.packages("testthat") # required only for tests
source("run_tests.R")
source("tests/core_checks.R")
```

`tests/test-bedtime-intervals.R` covers local date-time parsing, date/time column joining, interval roles and clipping, overlap and order errors, circular-mean baselines across midnight, column proposals, the new defaults, and the full Shiny workflow with an Actiware-style REST/SLEEP export, including reproduction of the exported analysis. Logs from this release are in `verification/`.

---

# Earlier release notes (1.3.4)

## English KSS labels on the composite graph

The composite index graph now always shows all nine English KSS descriptions on the right axis, including when observation overlays are off. The graph is taller and long labels wrap onto two lines. Horizontal zoom remains available; the vertical scales stay fixed so the descriptors remain aligned.

The labels use the [English nine-level KSS anchors](https://pmc.ncbi.nlm.nih.gov/articles/PMC5511283/). Their positions correspond to `(KSS - 1) / 8`; this is a reference scale and does not validate a conversion of the weighted composite to KSS. The model and index calculations are unchanged.

## Fix for names such as KSS_10 and KSS_10_SD

Under the empty measurement-times box, click **Use column numbers as clock hours**. This is an explicit confirmation that the numeric part represents an hour of day, not an occasion index. It applies to KSS, VAS and PVT. A proposal is shown before you click; existing nonempty/manual values are never replaced.

For the headers shown in the earlier screenshot (`KSS_10` through `KSS_24`, then `KSS_2` through `KSS_8`, then `KSS_10_SD` through `KSS_20_SD`), the recording hours become `10, 12, …, 44`. At an 08:00 study start, the elapsed hours are `2, 4, …, 36`.

`24` is interpreted as midnight (`0`) before unwrapping. Text such as `_SD` is not interpreted as a day marker or excluded automatically: the user must select the appropriate measurement columns and check their recording order. If `_SD` means a statistic rather than a measurement, exclude those columns.

Incomplete or invalid timing produces local validation feedback without the reported errors in the composite and cross-modality panels. Simulation and export remain blocked while an included dataset has invalid times; datasets are not silently omitted from analysis.

## Clock times for KSS, VAS and PVT

Explicit clock labels populate the compact measurement-times box automatically, in **file-column order**. A green note asks you to check the proposed values.

| Names | Proposed hours |
| --- | --- |
| `1h00`, `2h00`, `18h00` | `1, 2, 18` |
| `08:00`, `10:30`, `23:00`, `00:00` | `8, 10.5, 23, 0` |
| `KSS_9u_dag1`, `KSS_21u_dag1`, `KSS_3u_dag2` | `9, 21, 3` |
| `11PM`, `12AM`, `1:30am`, `12PM` | `23, 0, 1.5, 12` |
| `Base7h30`, `Base9h` | `7.5, 9` |
| `T1`, `T2`, `T3` or `VAS_10`, `VAS_12` | No proposal |

If any selected name lacks an unambiguous explicit clock, enter the values in the box. Use decimal hours (`10.5` for 10:30). There must be exactly one finite, nonnegative value per selected measurement column. There is no automatic spacing fallback. Bare numeric suffixes require the explicit one-click confirmation described above; `T1`, `T2` remain uninterpreted.

Typed values and existing/restored time-input values are preserved. Changing study start or requesting another proposal does not overwrite nonempty text. A new dataset permits new proposals; changed column selections refresh untouched proposals, but preserve manual entries for you to correct.

A strict decrease adds 24 hours to that and every later value: `22, 0, 2` becomes `22, 24, 26`. A notification lists the adjusted recording hours. Equal consecutive times are rejected. For known successive days at the same hour, enter explicit hours such as `8, 32`.

Elapsed observation time equals unwrapped recording hours minus the study-start clock. At study start 08:00, `22, 0, 2` maps to elapsed hours `14, 16, 18`. Plots, comparisons, automatic duration suggestions and exports use these mapped times.

The first measurement's day is not inferred: a first time before study start is shown as outside coverage. If it belongs to the next day, enter its hour beyond 24. Day labels alone do not resolve equal times.

## Other corrections and scope

Study-start validation is caught and shown within the app while incomplete or invalid datetimes are being edited. Date-based mode still requires an explicit UTC offset.

Standalone exports now include the agreement-statistic helper and null-default operator used by this file's calculations. This fixes missing dependencies when reproducing an analysis in a fresh R session.

The supplied file's numerical, comparison and calibration functions are preserved, apart from those export changes and replacement of obsolete time-parsing helpers. Existing unvalidated/custom-parameter and composite-index labels remain. This release verifies the requested software changes; it is not a new scientific validation of the supplied file's added calibration or agreement methods.

### Dataset access limitation

macOS denied access to the referenced CSV even after read permission was granted. Regression tests therefore use synthetic data with the KSS headers visible in the earlier screenshot, not the contents of that CSV. The actual file can be checked separately if it is attached in the chat or copied to an accessible workspace folder.

