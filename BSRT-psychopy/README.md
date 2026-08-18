# BSRT (PsychoPy) — Behavioral Sleep Resistance Task

The same OSLER-paradigm task as `../BSRT` (browser) and `../BSRT-desktop`
(Electron), as a PsychoPy experiment. Same paradigm, same scoring, same CSV
columns — so data from all three builds `rbind`s and the same analysis reads it.

**Read this before running participants: [what is verified, and what is
not](#what-is-verified-and-what-is-not).** An earlier revision was confirmed
running on a real PsychoPy installation; the restyled screens added since have
been exercised only against a stub.

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

### Both paradigms

The first dialog chooses between them, because they need different defaults and
a modal dialog cannot re-tick a box when the mode changes the way the browser
build does.

| | **BSRT / OSLER** | **PVT** |
|---|---|---|
| Stimulus | a red dot | a counting millisecond timer |
| Interval | fixed, 3000 ms | varies: 2, 4, 6, 8, 10 s |
| Sleep-onset criterion | **on** (7 consecutive misses) | **off** |
| Ends | at the criterion, or the ceiling | at the ceiling |

The PVT schedule is a fresh random **permutation** of the interval set per 30 s
block — random in order, balanced in composition, so no block comes out
all-short or all-long. The interval set is editable in the dialog
(`isi_set_s`, e.g. `2/4/6/8/10`); a set that does not sum to the block length
cannot be balanced, and the generator says so by recording
`sampled_with_replacement` in `schedule_method` rather than pretending.

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

## The look

The browser build is styled with CSS; PsychoPy has none, so the same design is
rebuilt in `bsrt_ui.py` out of rectangles and text — the same palette taken
from `BSRT/styles.css` (`--bg`, `--panel`, `--accent`, the three band colours),
the same titled-card layout on a dark ground, muted secondary text, and a key
hint along the bottom.

What that changes in practice:

- **The KSS is a scale**, not nine numbered lines: a row of boxes, the ends
  always labelled so the direction is unambiguous, the selected anchor spelled
  out underneath. Answer with 1–9, or arrow keys and Enter.
- **The countdown** is a number in an accent ring rather than a bare digit.
- **Sleep onset** takes over the screen in the danger colour, with the time it
  happened, and holds until the alarm is silenced.
- **The results screen** lists the headline numbers and, when a comparison was
  available, the norm bands — each carrying a colour *and* a coloured edge *and*
  the word "worse" or "better", exactly as the browser table does, so it
  survives greyscale and red/green colour blindness.

Stimuli are created once and reused rather than rebuilt per frame, which is the
usual reason a PsychoPy screen stutters. Nothing in `bsrt_ui.py` decides
anything: the task logic is in `bsrt_task.py` and the scoring in
`bsrt_core.py`, so a change to the look cannot change a number.

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

**The files are written the moment the trial ends**, before the sleep-onset
banner, the closing KSS or the results screen run. Those are screens, and a
screen can fail; none of them is worth a participant's data. If one does fail,
the traceback is printed along with the paths that were already saved, and the
run continues instead of ending on something that looks like a lost session.
(The closing KSS answer is added by rewriting the same four filenames, so
nothing accumulates.)

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

The task screens run in **English, French, Dutch and German**, set by the
`language` field in the first dialog. Every string on them — instructions, the
KSS and its anchors, the countdown, the sleep-onset banner and the results —
comes from `strings.json`, which `tools/gen_strings.py` extracts from
`BSRT/i18n.js`. All three builds therefore show the same wording, and
`tests/test_psychopy_path.py` runs a whole trial in French and reads back every
piece of text that reached the screen to prove nothing English leaks through.

**The settings dialogs stay in English.** They are PsychoPy's own
`DlgFromDict`, which labels each row with the setting's name — the same names
used in this README, in `DEFAULTS`, and in the exported columns. They are the
experimenter's screens, filled in before the participant sits down, and the
participant never sees them. Translating them would mean translating field
identifiers, which would then no longer match the documentation. Say the word
if you would rather have them translated anyway; it is a contained change.

**The French, Dutch and German KSS anchors have not been through formal
translation verification.** They belong to a validated instrument — check them
against an officially validated translation for your language before
publishing, and correct `i18n.js` (not `strings.json`, which is generated).

## What is verified, and what is not

Four test suites, run by `python3 tests/run_all.py`:

| Suite | What it establishes |
|---|---|
| `test_equivalence.py` | the Python scoring and `scoring.js` agree to 1e-9 across random trials, edge cases, five schedules and 120 RNG draws |
| `test_task.py` | the trial loop: stimulus counts, exact RTs, onsets within one frame with no drift over 60 epochs, the sleep-onset stop, dropped-frame detection, PVT vs BSRT, aborts |
| `test_headers.py` | all four CSVs match `BSRT/app.js` column for column |
| `test_export.py` | a whole simulated trial writes four files with every row the width of its header, and the right values in them |
| `test_psychopy_path.py` | the PsychoPy branch — window, keyboard, KSS, countdown, alarm-and-silence, abort, export — runs end to end |

**What is NOT established: none of this runs against a real PsychoPy.** PsychoPy
cannot be installed in the environment where this is written (its `esprima`
dependency fails to build there), so `test_psychopy_path.py` runs it against a
stub in `tests/stub_psychopy/`. That stub catches typos, wrong attribute names
and flow bugs in *this* code, but it is written from the same reading of the
PsychoPy documentation as the code it tests — so a misread API would be
reproduced in both and pass.

An earlier revision **was** reported working on a real installation, which
settles the window, keyboard and timing calls. The screens in `bsrt_ui.py` came
after that and have not been seen on real hardware: they use `Rect`,
`TextStim` positioning and `colorSpace='rgb255'`, none of it exotic, but
unconfirmed.

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
