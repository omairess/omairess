# DANCE: working rules

## The exported R script must replicate every analysis without the app

Students must be able to reproduce every result from the exported code alone
(Data Export → R script, `generate_analysis_code()` in `server/90_export.R`).
So **every change that can alter a number, a table or a test must reach the
script in the same commit**. That covers:
- a new option or setting, and any changed default;
- a new estimator, outcome or comparison;
- a new module.

Routine for every such change:

1. **Put it in the script.** Pass the new setting to the call the script makes,
   or add the module's section.
2. **Call the app's own functions.** Put the estimator in a pure helper in
   `server/08*_helpers_*.R` and write it into the script with `emit_kernel()`,
   which also emits every helper it calls. Never retype an estimator in the
   script.
3. **Write constants as values.** Settings, data vectors and options go into
   the script as literal values (`dance_export_vec()`, `deparse()`), never as
   references to app state.
4. **Prove it with a round-trip test.** Run the app, generate the script, run
   the script in a fresh environment on the same data (`dance_input_data`), and
   check that its numbers equal the app's. Examples: `tests/cosinor_export_test.R`
   and `tests/dlmo_app_test.R`. A new option gets its own case.
5. **State what was done.** The script's comments carry the settings text, the
   same words as the summary and the report.

A change that only affects display (a colour, a click, a layout) needs no
script change; say so when reporting it.

## Other standing rules

- **No new procedure without validation.** Reuse an established package
  where one exists, and check for one first. When a calculation has to be
  done in-house, log it in `CALCULATIONS_LOG.md`: the rule, its source, its
  tests.
- **Nothing assumed silently.** Every default and every decision a run makes
  is stated in the summary, the report and the script.
- **Keep it simple.** Keep the output decluttered.
- **After each change:**
  1. run the tests;
  2. commit and push;
  3. build the zip for testing:
     `git archive --format=zip --prefix=DANCE/ HEAD -o DANCE.zip`, then check
     it with `unzip -tq`.
- **Running the tests:**
  - App tests run from `DANCE/` with `LANG=C.UTF-8 Rscript tests/<x>.R`.
  - Unit tests run with `testthat::test_dir("tests/testthat")`.
