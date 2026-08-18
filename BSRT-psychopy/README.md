# BSRT (PsychoPy) — Behavioral Sleep Resistance Task

The same OSLER-paradigm task as `../BSRT` (browser) and `../BSRT-desktop`
(Electron), as a PsychoPy experiment. Same paradigm, same scoring, same CSV
columns — so data from all three builds `rbind`s and the same analysis reads it.

**Read this before running participants: [what is verified, and what is
not](#what-is-verified-and-what-is-not).** This build has never been executed
against a real PsychoPy installation.

---

## Running it

```bash
python bsrt_psychopy.py
```

Or open it in the PsychoPy **Coder** view and press Run. A dialog collects the
participant, the task settings and where to save; the trial then runs
fullscreen.

There is no Builder `.psyexp`. The task is a timed loop with a stopping rule
that Builder expresses awkwardly, and it has to stay identical to the other two
builds — a Builder version would be a fourth implementation to keep in step.

### Controls

| Key | |
|---|---|
| space, or any other key | respond |
| escape | abort the trial — **the data collected so far is still saved** |
| s | silence the sleep-onset alarm |

`escape` and `s` are experimenter keys and are never recorded as responses. If
you need `s` to be a response key, change `SILENCE_KEYS` in `bsrt_psychopy.py`.

### Without a screen

```bash
python bsrt_psychopy.py --simulate --minutes 3            # BSRT
python bsrt_psychopy.py --simulate --minutes 3 --mode pvt # PVT
```

Runs the whole pipeline — schedule, trial loop, scoring, all four CSVs — on a
virtual clock with scripted responses. Useful for checking an analysis script
against a file of the right shape before anyone is tested, and it is how this
build is exercised in CI.

---

## Why a PsychoPy build

Timing. This build takes a stimulus onset from the **vsync timestamp of the
frame that actually carried it**, and a reaction time from the **keyboard's own
hardware timestamp** against that onset. The browser build can only approximate
both, and says so at length in its own README.

| | browser | desktop (Electron) | PsychoPy |
|---|---|---|---|
| Onset | `requestAnimationFrame` estimate | frame clock | vsync timestamp |
| Response | DOM event stamp, corrected | DOM event stamp, corrected | keyboard hardware stamp |
| Clock resolution | 0.1–1 ms, coarsened | 0.1–1 ms | sub-millisecond |
| Install needed | none | one download | PsychoPy |

If PsychoPy is already set up in your lab, this is the most accurate of the
three. If it is not, the desktop build is close and needs no Python.

### Epoch boundaries come from the schedule, not from the last onset

Each stimulus is shown on the first frame at or after its **intended** onset,
and each epoch closes at its intended boundary. Closing an epoch a fixed
interval after the frame it happened to appear on lets each epoch's presentation
error feed into the next, and it compounds — the first version of this loop
drifted 25 ms per epoch and put the sixtieth stimulus of a 3-minute test 1.5 s
late. The measured onset of every stimulus is recorded next to its intended one,
so presentation error is a number you can look at rather than an assumption.

---

## What is shared with the other builds

Nothing that matters is retyped. Four things are generated or checked against
the browser build:

| Here | Comes from | Checked by |
|---|---|---|
| `bsrt_core.py` | a port of `BSRT/scoring.js` | `tests/test_equivalence.py` |
| `norms.json` | `tools/gen_norms2.py`, same run that writes `norms.js` | byte comparison |
| `strings.json` | `tools/gen_strings.py`, extracted from `BSRT/i18n.js` | regenerate |
| CSV columns | `BSRT/app.js` | `tests/test_headers.py` |

`bsrt_core.py` is deliberately a line-for-line port rather than an independent
implementation: a second reading of the spec would drift, and the drift would
show up as a systematic difference between data collected in different rooms.
The equivalence test drives the Python and the JavaScript with the same inputs
— including the seeded schedule generator, so the same seed produces the same
schedule in all three builds — and asserts they agree to 1e-9.

To regenerate after changing the shared sources:

```bash
BSRT_NORM_DIR=/path/to/workbooks python3 tools/gen_norms2.py --write
python3 tools/gen_strings.py
python3 BSRT-psychopy/tests/run_all.py
```

---

## Output

The same four CSVs as the other builds, written to the output directory:

```
bsrt_<id>_<session>_t<n>_raw.csv          one row per stimulus — the irreplaceable one
bsrt_<id>_<session>_t<n>_perminute.csv    one row per minute
bsrt_<id>_<session>_t<n>_summary.csv      one row for the trial
bsrt_<id>_<session>_t<n>_norms.csv        one row per normative variable
```

Columns are identical to the browser build's, in the same order. Columns that
describe a browser (`device_browser`, `screen_w`, `page_host`,
`cross_origin_isolated` and the rest) are written **empty** rather than dropped:
the shape stays identical so the files still `rbind`, and empty means "this
build cannot know", which is not the same as zero.

Columns this build fills differently:

| Column | Here |
|---|---|
| `device_platform` | `psychopy` |
| `rt_source_setting`, `rt_source_used` | `psychopy_flip_to_kb` |
| `rt_event_ms`, `rt_handler_ms`, `input_delay_ms` | **empty** — there is one clock, so there is no second reading to correct against |
| `refresh_hz_measured`, `frame_interval_ms` | measured at startup |
| `onset_quantisation_ms` | the frame interval: on a frame-locked build that *is* the quantisation |
| `frames_trial`, `dropped_frames_trial` | counted during the trial |

The schedule seed is exported as `schedule_seed`, so the exact sequence a
participant saw can be regenerated from the CSV alone.

---

## Norms

Identical to the other builds: 435 sessions from two pooled studies, binned by
hour of day and cumulative test length, with a comparison written to the norms
CSV. **Read `../BSRT/README.md` on what the reference is and what it is not** —
two samples of unequal sleepiness, a survivor bias at the longer lengths, and
hour bins holding about a dozen sessions each.

This build writes the comparison to the CSV but does not draw the coloured
table on screen; the experimenter sees a short end-of-trial summary instead.

---

## Languages

The on-screen text and the KSS anchors are available in English, French, Dutch
and German, set with the `language` field in the dialog. They are extracted from
`BSRT/i18n.js`, so all builds show the same wording.

**The French, Dutch and German KSS anchors have not been through formal
translation verification.** They belong to a validated instrument — check them
against an officially validated translation for your language before publishing,
and correct `i18n.js` (not `strings.json`, which is generated).

---

## What is verified, and what is not

Four test suites, run by `python3 tests/run_all.py`:

| Suite | What it establishes |
|---|---|
| `test_equivalence.py` | the Python scoring and `scoring.js` agree to 1e-9 across random trials, edge cases, five schedules and 120 RNG draws |
| `test_task.py` | the trial loop: stimulus counts, exact RTs, onsets within one frame with no drift over 60 epochs, the sleep-onset stop, dropped-frame detection, PVT vs BSRT, aborts |
| `test_headers.py` | all four CSVs match `BSRT/app.js` column for column |
| `test_export.py` | a whole simulated trial writes four files with every row the width of its header, and the right values in them |
| `test_psychopy_path.py` | the PsychoPy branch — window, keyboard, KSS, countdown, alarm-and-silence, abort, export — runs end to end |

**What is NOT established: this build has never been run against a real PsychoPy
installation.** PsychoPy could not be installed in the environment where it was
written (its `esprima` dependency fails to build there), so
`test_psychopy_path.py` runs it against a stub in `tests/stub_psychopy/`. That
stub catches typos, wrong attribute names and flow bugs in *this* code, but it
is written from the same reading of the PsychoPy documentation as the code it
tests — so a misread API would be reproduced in both and pass.

Before using this for real data, run it once on a real machine and check:

1. **It starts at all.** Import errors and `visual.Window` arguments are the
   likeliest first failure.
2. **`win.flip()` returns a timestamp.** The engine treats the return value as
   the time the frame appeared. If your PsychoPy version returns `None`, the run
   will fail immediately and loudly rather than silently mis-time anything.
3. **Key timestamps share a clock with the flips.** Press the key in time with a
   visible stimulus and check the recorded RTs are plausible — a constant offset
   of hundreds of milliseconds, or negative RTs, means `kb.tDown` and `flip()`
   are on different clocks. Install the `psychtoolbox` backend if it is missing;
   PsychoPy falls back to a less accurate keyboard without it.
4. **`refresh_hz_measured` in the summary matches your monitor**, and
   `dropped_frames_trial` is near zero on a run where nothing else was
   happening.

Report anything that fails and it can be fixed against real behaviour rather
than against my reading of the docs.

---

## Requirements

PsychoPy 2021.2 or newer, which brings its own Python. `psychtoolbox` is
strongly recommended — it is what gives the keyboard its hardware timestamps.
The scoring, scheduling and export code is **stdlib only**: no numpy, no pandas,
nothing to install beyond PsychoPy itself.
