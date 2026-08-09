# Rebuilding the normative tables

`BSRT/norms.js` and `BSRT-desktop/renderer/norms.js` are generated, not written
by hand. This directory holds what generates them.

```bash
pip install openpyxl
BSRT_NORM_DIR=/path/to/workbooks python3 tools/gen_norms2.py           # dry run
BSRT_NORM_DIR=/path/to/workbooks python3 tools/gen_norms2.py --write   # emit both copies
```

The dry run prints the bin sizes, the censoring table, and a side-by-side of the
two studies, and touches nothing. `--write` writes both copies of `norms.js`;
run `npm run check:shared` in `BSRT-desktop` afterwards to confirm they are
byte-identical.

## The source workbooks are not in this repository

Two files are expected in `BSRT_NORM_DIR`:

| File | Study |
|---|---|
| `BSRT_CTRL_GROUP.xlsx` | control group, 291 sessions, 8-minute protocol |
| `BSRT_outcome_calculations_IH_2024.xlsx` | overnight study, 144 sessions, 40-minute protocol |

They hold participant-level trial data including **dates of birth**, so they are
deliberately kept out of version control. Keep them under your data-protection
plan and point the environment variable at them.

## What the code does

`normlib.py` is a port of `normativeSummary()` from `scoring.js` — deliberately
so, because the reference values and a participant's own values have to be
computed by the same rule or every z-score carries a systematic offset. If you
change the scoring conventions in one, change them in the other and regenerate.

`gen_norms2.py` reads raw trials only. It never reads the workbooks' own
precomputed `Mean_*` / `SD_*` columns, which use a different convention (they
keep no-response timeouts in the RT average; the app excludes them). Those
columns were used once, as an independent check that the raw parsing was right,
and that check is recorded in the header of the generated file.
