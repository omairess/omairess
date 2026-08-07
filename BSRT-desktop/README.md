# BSRT (desktop) — Behavioral Sleep Resistance Task

An OSLER-paradigm sleep resistance task packaged as a native application for
**macOS and Windows**, with frame-locked stimulus presentation and a calibration
routine that measures and reports the timing characteristics of the machine it
is running on.

This is the accuracy-focused sibling of the browser version in `../BSRT`. The
paradigm, scoring, and CSV shape are the same. What differs is how time is
measured.

---

## Getting the app

### Easiest: download a build (no tools needed)

The repository builds ready-to-run apps on real macOS and Windows machines, so
nobody has to install anything to get one.

1. Open the repository on GitHub → **Actions** tab
2. Choose **Build BSRT desktop apps** → **Run workflow**
3. Wait about five minutes
4. Download **BSRT-macOS** (a `.dmg`) or **BSRT-Windows** (an installer and a
   portable `.exe`) from the artifacts at the bottom of the run

Pushing a version tag (`git tag v1.0.0 && git push --tags`) also publishes them
to a **Releases** page, which gives a permanent link you can send to
collaborators — they just click and download.

### Building it yourself

Only needed if you want to change the app. Requires
[Node.js](https://nodejs.org) 18 or newer.

```bash
cd BSRT-desktop
npm install
npm start          # run it without packaging
npm run dist:mac   # on a Mac  -> dist/*.dmg
npm run dist:win   # on Windows -> dist/*.exe
```

**macOS apps must be built on a Mac and Windows apps on Windows** — packaging is
not portable, which is exactly why the CI workflow exists.

### The unsigned-app warnings

These builds are not code-signed, because signing needs a paid Apple Developer
account and a Windows code-signing certificate. On first launch:

- **macOS:** right-click the app → *Open* → *Open*
- **Windows:** *More info* → *Run anyway*

Once per computer. If you are distributing to participants rather than running
the task yourself, budget for certificates — talking someone through a security
warning by phone is a poor start to a testing session. Set `CSC_LINK` /
`CSC_KEY_PASSWORD` (macOS) or `WIN_CSC_LINK` / `WIN_CSC_KEY_PASSWORD` (Windows)
as repository secrets and the CI build will sign automatically.

---

## Languages

The interface runs in **English, French, Dutch and German**, switchable at any
time from the selector in the top-right corner — including on the results
screen, which re-renders in the new language. The language a participant saw is exported as `language`
with every trial.

**Translation status:** English is the source. The other three were produced for
this app and have **not** been through formal translation verification. The
Karolinska Sleepiness Scale anchors in particular belong to a validated
instrument — before publishing, check them against an officially validated
translation for your language and correct `i18n.js` if they differ. Any string
missing from a translation falls back to English rather than appearing blank.

---

## Karolinska Sleepiness Scale

Optional, and independently schedulable: **not at all, before the task, after
the task, or both**. The standard 9-point form is used, with an anchor on every
step, presented one screen at a time between the setup screen and the countdown
(and again after the task, before the results).

Answers export as `kss_before` and `kss_after`, with `kss_when` recording what
was requested. A trial with the KSS off exports empty cells rather than zeros,
so "not asked" is never confused with "answered 0".

---

## False starts and response integrity

A **false start** is a response between 0 and 100 ms after stimulus onset — too
fast to be a genuine reaction. It applies to both paradigms, is flagged per
trial as `false_start`, totalled as `total_false_starts`, and is what the
"remove false starts" correction removes. The threshold is configurable and
exported as `false_start_threshold_ms`.

Continuous or repeated pressing is a way to fake alertness: hold or tap
constantly and every stimulus gets a fast "response". The app measures three
independent signals and reports all of them, so you can see *why* a trial was
flagged:

| Column | Meaning |
|---|---|
| `total_presses` | every keypress in the trial, including presses outside any stimulus |
| `extra_presses` | presses beyond the first in an epoch |
| `burst_max` | the most presses inside any 1-second window |
| `rapid_pairs` | consecutive presses closer together than 200 ms |
| `cheating_suspected` | 1 if any signal exceeded its threshold |
| `cheating_reasons` | plain text naming which, and by how much |

**This flags a trial for human review. It never alters scoring and never
excludes data by itself** — a startled double-tap is not cheating, and the
decision is yours. Thresholds are exported with the data.

---

## Alarm

When the sleep-onset criterion is reached, a **continuous pulsing alert** starts
and **keeps sounding until someone presses Stop**. A banner appears across the
top of the window with the stop button, on whatever screen is showing — so a
sleep-onset event cannot be missed by an experimenter who stepped away. It is on
by default, can be switched off, and the setting is exported as `alarm_enabled`.

The tone is synthesised (a square carrier gated by a 4 Hz square LFO) rather
than played from a sound file, so the app stays dependency-free and works
offline.

---

## The PVT stimulus

In **PVT mode the stimulus is a millisecond counter**, as on a PVT-192: it
appears and counts up from 0. Pressing stops it and the achieved reaction time
stays on screen as feedback before clearing. If no response comes, it stops at
the hit window.

The counter is **display only** — reaction time is still measured from stimulus
presentation to the keypress, never read off the digits. On the desktop build
the counter is redrawn on the frame clock, so it advances in exact frame steps
(0, 17, 33, 50 ms at 60 Hz) and redrawing cannot disturb epoch scheduling.

**BSRT mode keeps the red dot.** The two paradigms use their conventional
stimuli.

---

## Two modes

The task runs either paradigm; the only structural difference is how the
inter-stimulus interval is chosen.

| | **BSRT / OSLER** | **PVT** |
|---|---|---|
| Interval | fixed (3000 ms) | varies: 2, 4, 6, 8, 10 s |
| Ends | on the sleep-onset criterion, or at the ceiling | at the ceiling |
| Sleep-onset criterion | on | off by default |

Everything else — hit window, lapse threshold, scoring, corrections, per-minute
and total metrics, all three CSV exports — is identical, so the two are directly
comparable and the same analysis code reads both.

### How the PVT schedule is built

The specified intervals sum to exactly one block: 2 + 4 + 6 + 8 + 10 = **30 s**.
Each 30-second block is therefore a fresh random **permutation** of the five
intervals — random in order, balanced in composition. Every block delivers one
of each wait, so no block can come out all-short or all-long, and the
interval distribution is identical in every block of the test.

This is an interpretation of "vary randomly within each 30 s block". If you
intended independent sampling with replacement instead, say so — it is a
one-line change. The rule actually used is exported as `schedule_method`
(`block_permutation` or `sampled_with_replacement`), so it is never ambiguous
after the fact.

A custom interval set that does *not* sum to the block length cannot be
balanced. The generator says so by falling back to sampling with replacement and
recording that in `schedule_method`, rather than silently pretending to balance.

### The wait comes before the stimulus

Each interval **precedes** its stimulus, so the screen stays blank for one full
interval after the countdown rather than flashing the instant it ends. In BSRT
that is 3000 ms; in PVT it is 2, 4, 6, 8 or 10 s, drawn like any other interval,
so a participant cannot learn when the first one is coming.

A schedule of N intervals yields **N stimuli**: a 3-minute, 3000 ms BSRT
presents 60, and a 3-minute PVT presents 30, matching what the settings imply.
The last stimulus gets a response window as long as the longest interval in the
configured set — 3000 ms for BSRT, 10 s for PVT — the same as every other
stimulus, rather than being cut off the instant it appears. The trial therefore
runs a little past the requested ceiling: a 40-minute BSRT lasts 40 minutes and
3 seconds, and a PVT with the default set lasts its ceiling plus 10 s. That
overshoot is fixed by the settings alone, not by the random schedule, so every
participant still gets exactly the same duration and the same number of
stimuli whatever order the blocks were shuffled into.

Two columns describe each stimulus's timing: `isi_before_ms` is the wait that
preceded it (defined for the first stimulus too) and `epoch_isi_ms` is its
response window, the time until the next one.

A stimulus belongs to the block its **preceding interval** started in, so each
30-second block still holds exactly one stimulus per interval in the set.

### Reproducible schedules

Schedules come from a seeded generator, and the seed is exported as
`schedule_seed`. Leave the field blank for a random seed per trial, or type one
to reproduce a schedule exactly — useful for a fixed test order across
participants, or for reconstructing what a participant actually saw.

### What is recorded per stimulus

The raw CSV carries `block`, `epoch_isi_ms` (this stimulus's interval) and
`isi_before_ms` (the wait that preceded it). Reaction time depends on the
preceding interval in a PVT, so that column is what you need to model it.

### A consequence worth knowing

The response window is the whole epoch, so with a 10 s interval a very late
press is recorded against that stimulus with a reaction time of several
seconds. It is still scored a **miss** (`late_response = 1`) and never counts as
a hit, but if you want presses during the dead time treated as false starts
against the *next* stimulus instead, that is a different rule and would need
adding deliberately.

---

## Trial procedure

One epoch, in order:

1. A new stimulus starts every **3000 ms**, measured from the previous
   stimulus's onset.
2. The light stays on for up to **1000 ms** (the *hit window*).
3. **If the participant responds, the light goes out immediately.**
4. If they do not respond, it goes out by itself at 1000 ms.
5. The next stimulus still starts at exactly +3000 ms from this one's onset —
   **never from the response**. The epoch clock is independent of behaviour, so
   a trial lasts the same for a fast responder and a slow one.

### How each epoch is scored

| RT | Outcome |
|---|---|
| `RT <= 1000 ms` | **HIT** |
| `RT > lapse threshold` (and still a hit) | **LAPSE** — a slow hit |
| `1000 ms < RT < 3000 ms` | **MISS** — the raw RT is still recorded |
| no response | **MISS** |

A response slower than the hit window is a miss even though a key was pressed:
it does **not** reset the consecutive-miss run, so it counts toward sleep onset
exactly like a silent epoch. Its reaction time is still stored and exported
(`late_response = 1`), because raw data is never discarded.

A lapse is a *kind of hit*, counted inside the hit total, and never contributes
to the sleep-onset criterion.

**Lapse threshold** defaults to **500 ms**. For short protocols (under about
5 minutes) **355 ms** is conventional; the setup screen says which applies once
you set the test duration. The value used is exported as `lapse_threshold_ms`
so the definition travels with the data.

---

## Outputs

Three CSV files per trial, all long-format with participant, date, time and
trial keys on every row so files `rbind` cleanly in R.

**1. Raw — one row per stimulus.** Always exported, never filtered.
`epoch_index`, `block`, `minute`, `onset_ms`, `epoch_isi_ms`, `isi_before_ms`,
`responded`, `rt_ms`, `rs_per_sec`, `outcome`, `lapse`, `late_response`,
`false_start`, `extra_responses`.

`epoch_index`, `block` and `minute` all start at **1**, not 0.

`block` and `minute` are taken from the planned schedule, not re-derived from
the measured onset: a stimulus is assigned to the block and minute in which its
*preceding interval* started. That matters at the boundaries. Bucketing on the
onset instead puts a stimulus landing exactly on the minute mark into the
following minute, and since onsets span (0, ceiling] rather than [0, ceiling),
a 2-minute test came out as 19 / 20 / 1 across *three* minutes — the third
holding a single trial whose mean, velocity and acceleration meant nothing. It
also has to come from the schedule rather than the measurement, because
presentation jitter of a millisecond (or a whole frame, on the desktop build)
is otherwise enough to tip a boundary stimulus into the wrong bucket. On the
scheduled rule a 2-minute, 3000 ms test reads exactly 20 / 20.

**2. Per-minute — one row per elapsed minute.**
`trials`, `hits`, `misses`, `lapses`, `hit_ratio`; average / median / SD /
10% fastest / 10% slowest for both **RT** and **RS**, each in raw and corrected
form; plus `velocity_rs`, `acceleration_rs`, `velocity_rt`, `acceleration_rt`.

**3. Summary — one row per trial.**
Participant information, `sleep_onset_ms`, `total_trialrun`, `hit_ratio`,
`total_miss`, `total_lapse`, the error profiles `ep_1_2` / `ep_3_6` /
`ep_7plus`, every whole-test RT and RS metric (raw and corrected), the
regression slopes, and the full protocol configuration.

### Definitions used

**Reaction speed** is `RS = 1000 / RT`. It is undefined for `RT <= 0`; such
trials are excluded from RS and counted in `rs_undefined` rather than silently
dropped.

**10% fastest / slowest** are means of the extreme decile (at least one trial).
For RT the fastest responses are the *low* tail; for RS they are the *high*
tail. The two are computed independently from their own distributions.

**Error profiles** count **runs** of consecutive misses, binned by run length —
EP1–2, EP3–6, EP7+ — not the number of missed trials. A run of 7 ends the test,
so `ep_7plus` is normally 0 or 1.

**Corrections** are configurable: none, anticipations only (`RT < 100 ms`),
outliers only (beyond 2 SD from the mean), or both. Anticipations are removed
**first**, so they cannot inflate the mean and SD that define the outlier
bounds. Bounds are computed within whichever block is being summarised — a
per-minute correction uses that minute's own mean and SD, a whole-test
correction uses the whole test's. Raw columns are never corrected.

**Velocity and acceleration** are the first and second differences of average
reaction speed across successive minutes. Velocity is undefined for minute 1 and
acceleration for minutes 1–2; these are reported as empty, not zero. Positive
velocity means the participant is speeding up. `rs_slope_per_min` is the
least-squares slope of reaction speed against minute across the whole test.

### Participant identifiers

Birth dates are entered from three dropdowns — day, month, year — rather than
typed. Free text produced `1985-03-07`, `07/03/1985` and `March 1985` from
different testers on the same study, and the resulting column could not be
parsed. There is nothing to type, so an out-of-format answer cannot be
produced; the stored and exported value is always ISO `YYYY-MM-DD`. Month names
follow the interface language. The day list follows the month, including leap
years, so 31 February is never on offer. A half-finished date (say a month and
a year but no day) is refused with an explanation rather than exported as a
fragment; leaving all three blank stays a valid answer for an anonymised
participant.

Dropdowns are used in preference to a native calendar popup because a birth
date sits decades back: picking 1985 from a list is one click, where a calendar
means paging through hundreds of months.

### Returning participants

Everyone who completes a trial is added to a small local roster, and a second
visit starts by choosing them from the **Returning participant** dropdown —
their ID, name, birth date, address and educational level are filled in, and
the trial number advances to the next one free for that participant and
session. Typing a known ID has the same effect as you type. Nothing has to be
entered twice, which is where a `P001`/`POO1` or a transposed birth year
normally creeps in and quietly splits one participant into two.

The roster is kept separately from the trial data, so clearing stored trials
does not lose it, and an installation that already holds trials has its roster
rebuilt from them the first time it starts.

### Protection against accidental overwriting

Recall alone would make overwriting easier, so two checks sit in front of the
Start button. Neither can be dismissed by clicking past it: the trial does not
begin until the conflict is resolved.

**Changed identity.** If any identity field no longer matches the saved record,
the difference is listed field by field — old value struck through, new value
beside it. Clearing a stored value counts, because losing a birth date by
accident is as damaging as replacing it with a wrong one. Two ways out:
*Restore saved details* puts the record back, or a tick box confirms that the
new details are the correct ones and the record should be updated. The tick
applies to one decision only — change the details again and it resets, so a
confirmation can never carry over to a conflict nobody looked at.

**Reused trial number.** If this participant, session label and trial number
already have data, the existing trial is named with its date and time. *Use the
next free number* moves on with one click; a tick box allows a deliberate
repeat. A different session label starts its own numbering.

A participant record is created or updated only once a trial has actually been
recorded, so an abandoned setup screen never changes the roster. Details filled
in automatically are also withdrawn automatically if the ID is typed on past a
match, so a new participant cannot inherit someone else's name and birth date —
but anything entered by hand is never cleared.

### Privacy

Name, address and birth date are recorded and written into every export. They
are directly identifying. For research use, prefer a pseudonymous participant ID
and hold the identity key separately under your data-protection plan — the app
warns about this on the setup screen. Date and time of testing are captured
automatically.

### Scoring code

All derived metrics come from `scoring.js`, which is byte-identical between the
browser and desktop builds so the two cannot diverge.
`npm run check:shared` in `BSRT-desktop` asserts that `scoring.js`,
`i18n.js` and `participants.js` are identical in both builds, and the
packaging workflow refuses to build if they are not.

---

## Preliminary norms

Results are compared against a **preliminary** control reference of 291 test
sessions from 12 participants, binned by the **hour the session started** and by
**cumulative test length**. A 3-minute test is compared against the first 3
minutes of the normative protocol, not against its whole length — the reference
was recorded as a single 8-minute test and each length bin is a prefix of it.

Each of 21 variables gets its value, the reference mean ± SD for that hour and
length, and a signed distance in SDs where **positive always means worse**,
whichever direction is bad for that variable. Colouring follows that number:

| Band | Meaning |
|---|---|
| green | within 1 SD of the reference, or better than it |
| orange | more than **1 SD worse** |
| red | more than **2 SD worse** |

Colour is never the only signal — each band also carries a left border and the
word *worse* or *better*, so the table survives greyscale printing and red/green
colour blindness.

### When no comparison is offered

The reference sessions are BSRT runs at a 3000 ms interval, which is 20 stimuli
a minute. The panel refuses, and says which condition failed, when the trial is
a **PVT** (about ten stimuli a minute), when the **interval is not 3000 ms**,
when the **start hour** cannot be read, or when the test did not complete a
**full minute**. Trial, miss and lapse counts would otherwise be compared
against a reference measuring something else. A test longer than 8 minutes is
compared over its first 8 minutes only; there are no norms beyond that.

### Limits you should read before quoting a z-score

- **They are preliminary.** Each hour bin holds only **9–15 sessions**, and the
  same participant contributes to many bins. The observations are not
  independent, and every SD is itself estimated from about a dozen numbers.
- **The reference protocol was 8 minutes.** Comparing a 3-minute test against
  the first 3 minutes is the right window, but the reference participants knew
  they were settling in for 8 minutes. Someone expecting a short test may pace
  themselves differently, so **comparisons below 8 minutes can be biased for
  that reason alone**. The app flags this whenever the window is shorter.
- **Some reference cells have almost no spread.** Every control session scored
  20 hits out of 20 in the first minute, so a difference of one or two trials
  can produce a double-digit SD figure. Where the reference SD is exactly zero
  no z-score is computed at all; the value is marked as outside the reference
  range instead of dividing by zero.
- **The comparison uses the reference workbook's definitions, not the app's.**
  The workbook counts every response however slow and sizes a decile with
  `ceil(10%)`; the app treats a response slower than the hit window as a miss
  and uses `round(10%)`. Across the 291 reference sessions those two differences
  move the RT mean by a median of 4.8 ms, so the comparison re-scores the trial
  the workbook's way. **The numbers in the norms panel will not always match the
  tables above it**, which is deliberate — the export keeps both.

### Verification

The embedded table was checked against the source workbook two ways: all 4032
published cells (mean and SD) round-trip exactly through `norms.js`, and
recomputing every published mean from the raw trial values reproduces it, apart
from one session missing a single trial row that moves eight cells by at most
0.9 ms. The app's own re-scoring was then checked against an independent
implementation over 2160 fields with no mismatches.

### Export

A fourth CSV, **norms comparison**, is long-format with one row per variable:
`norm_value`, `norm_ref_mean`, `norm_ref_sd`, `z_worse`, `band`, plus the hour,
window, reference n, and the `norm_below_protocol` / `norm_truncated` flags. It
`rbind`s with the other exports.

---

## Timing: what this actually fixes

This is the part that matters, so it is stated without marketing.

### What the browser version could not do

- Stimuli were scheduled with `setTimeout`, which is best-effort and drifts
  under load.
- The stimulus onset time was recorded as "when the code set the CSS class",
  which is not when the pixel changed. The gap is up to a full refresh interval
  plus compositing, and it varies.
- Responses were timestamped when JavaScript got round to handling the event,
  which includes any queuing delay in the renderer.
- Nothing knew the display's refresh rate, so none of the above could even be
  quantified.

### What this version does

1. **Frame-locked scheduling.** A single `requestAnimationFrame` loop drives the
   whole trial and counts real vsync intervals. An epoch is a fixed number of
   frames, so presentation is periodic by construction and cannot drift relative
   to the screen. Verified: across a trial, successive stimulus onsets stayed
   within a fraction of one frame of the nominal interval, with no accumulation.

2. **Onset times from the frame clock.** The recorded onset is the timestamp of
   the frame the stimulus was painted into, not the moment a class was added.

3. **Latency from measured timestamps.** Sleep latency is read off the observed
   frame times, never recomputed as `epoch index × nominal ISI`. The calibrated
   frame interval carries a small estimation error, and multiplying it by an
   epoch index would let that error accumulate across a 40-minute trial.

4. **Responses timestamped by the OS.** `event.timeStamp` carries the time the
   platform delivered the event, on the same clock as `performance.now()`.
   Reading it instead of calling `performance.now()` in the handler skips the
   renderer's queuing delay. The app *measures* that delay during calibration,
   reports it, and verifies the assumption — if a future Chromium changed the
   timebase, the app detects it, falls back to handler time, and records that it
   did so rather than silently reporting wrong numbers.

5. **No throttling.** `backgroundThrottling: false` plus the background-timer
   and renderer-backgrounding switches. Without these, Chromium quietly throttles
   `requestAnimationFrame` when it decides a window is idle or occluded, which
   would corrupt a long trial with no visible sign.

6. **Display sleep blocked at the OS level** via `powerSaveBlocker`, which is far
   more reliable than the web Wake Lock API — a necessity in a task whose whole
   point is that the participant stops moving.

### Variable intervals and dropped frames

Stimuli fire on the first frame at or after their intended time, compared
against the real frame clock. Converting intended times into frame *indices*
using the calibrated interval would multiply any error in that estimate by the
elapsed time — a 16.70 ms estimate against a true 16.667 ms drifts 0.2%, about
five seconds across a 40-minute trial. Measured against the frame clock instead,
presentation still lands on a frame boundary, the error stays within half a
frame, and nothing accumulates: in testing, the final stimulus of a 100-stimulus
schedule was 2.4 ms from its intended time.


Each interval is rounded independently to a whole number of frames, and the
boundaries accumulate, so a variable schedule is delivered as precisely as a
fixed one — measured deviation was under 0.1 ms per interval in testing.

**Dropped frames cost one stimulus up to one frame interval, and nothing
after it.** Because each onset is scheduled against its own absolute intended
time rather than by counting frames forward, a frame lost at one boundary does
not push the next one late. The delay lands in that epoch's recorded
`achieved_isi_ms`, and never in the reaction time — RTs are measured from the
timestamp of the frame on which the stimulus actually appeared, not from when
it was meant to appear.

What matters is therefore the **rate**, not the count, and the results panel
now reports both (`dropped_frames_trial`, `frames_trial`, `dropped_rate_trial`
in the export). A bare count cannot be read: 97 dropped frames is 0.03% of a
40-minute trial at 120 Hz and 1.3% of a one-minute one. Under about 1% is
ordinary operating-system scheduling noise and needs no action.

Zero is not the expectation on a real machine, and a 120 Hz display is held to
a stricter standard than a 60 Hz one — each frame has 8.3 ms of budget instead
of 16.7, so a hitch invisible at 60 Hz is counted at 120. Drop counts are
therefore not comparable across machines with different refresh rates; the
rate is. On Macs with ProMotion, note that macOS varies the panel rate with
on-screen content, and a task that holds a static screen between stimuli is
exactly the case adaptive refresh is built to throttle.

To reduce them: run full-screen on a single display, quit other applications,
turn on Do Not Disturb, and keep the machine on mains power — macOS and Windows
both throttle harder on battery.

### What remains, and cannot be fixed in software

**No software on a general-purpose OS can observe photon emission or physical
key travel.** What is left after the above is:

| Source | Magnitude | Character |
|---|---|---|
| Display pipeline depth (frames buffered between paint and panel) | 1–2 frames (17–33 ms at 60 Hz) | near-constant per machine |
| Panel response / pixel transition | ~1–20 ms, panel-dependent | near-constant |
| Input hardware polling — wired USB | 1–8 ms | jitter |
| Input hardware — Bluetooth keyboard | 10–30 ms, variable | **jitter, avoid** |
| OS → application event delivery | typically <1–3 ms | measured and removed |

So: **jitter is reduced to roughly the input polling interval — single-digit
milliseconds on a wired USB keyboard — and a constant offset of very roughly
20–60 ms remains, unmeasurable in software.**

That constant is the honest limit. It is also the least damaging kind of error:
it is identical for every trial on a given machine, so it cancels in
within-participant and within-device comparisons. It does *not* cancel when
comparing raw reaction times across different hardware, which is why every
number needed to characterise it is exported with the data.

**Sleep latency — the primary outcome — is unaffected by all of this.** It is
scored in 3-second epochs; a few tens of milliseconds is four orders of
magnitude below the measurement grain.

### Measuring the constant, if you need it

Set *Presentation pipeline* and *Measured hardware offset* under **Timing &
calibration**, and enable the **photodiode patch** — a white square in the
top-left corner that flips with the stimulus.

1. Tape a photodiode over the patch (Black Box Toolkit, Cedrus StimTracker, or
   any oscilloscope-plus-photodiode rig).
2. Have the rig also register a keypress, ideally via a response box that closes
   the same circuit.
3. Run a few dozen stimuli. The measured interval between light onset and key
   closure, minus the reaction time the app reports, is your constant offset.
4. Enter it as *Measured hardware offset*. Raw and corrected reaction times are
   both exported; nothing is overwritten.

Until you do this, leave the offset at 0 and report reaction times as raw. The
app never silently applies a correction you did not measure.

### The app tells you when the assumption is wrong

If more than 10% of responses come in under 100 ms, the results screen says so
and names the likely cause. Genuine reactions are essentially never that fast,
so a crop of them means the presentation-pipeline setting overestimates this
display's lag and every reaction time is shifted. Responses under 100 ms are
counted as **anticipations** and excluded from the reaction-time statistics
(standard PVT practice); the raw values are still exported per epoch so you can
apply a different rule in R.

---

## Before you run participants

Frame-locked timing assumes a **fixed** refresh rate. Several common settings
break that assumption, and the calibration check will show it as jitter:

- **Turn off variable refresh rate** — G-Sync, FreeSync, VRR — in the GPU
  control panel.
- **macOS ProMotion displays vary between 24 and 120 Hz.** Set a fixed refresh
  rate in *System Settings → Displays* rather than *ProMotion*.
- **Run on mains power.** Laptops throttle the GPU on battery, which shows up as
  dropped frames.
- **One display only.** Mirroring or extending can pin the compositor to the
  slowest attached screen. The app warns if it detects more than one.
- **Close other applications**, particularly anything doing video or screen
  capture.
- **Use a wired USB keyboard**, not Bluetooth. This is the single largest
  remaining source of *jitter* and it is entirely under your control.

Run the display check once on each machine before testing and keep the numbers.
If it grades `poor`, fix the machine rather than the data.

---

## What gets recorded

Every trial carries its own timing provenance, so a reviewer can see exactly
what the numbers mean. The summary CSV includes, alongside the usual protocol
and performance columns:

`refresh_hz_measured`, `refresh_hz_reported`, `frame_interval_ms`,
`frame_jitter_mad_ms`, `onset_quantisation_ms`, `dropped_frames_calibration`,
`dropped_frames_trial`, `calibration_grade`, `input_time_source`,
`input_dispatch_median_ms`, `presentation_offset_frames`, `hardware_offset_ms`,
`isi_achieved_ms`, `n_rt_valid`, `anticipations`, `platform`, `os_release`,
`electron_version`, `chrome_version`, `display_scale_factor`.

Note `refresh_hz_measured` versus `refresh_hz_reported`: the first is what the
compositor actually achieved, the second is what the OS claims. When they
disagree, trust the measured one and investigate the machine.

`isi_achieved_ms` is the real inter-stimulus interval, which is the requested
interval rounded to a whole number of frames — at 60 Hz, a requested 3000 ms
becomes 2999.7 ms or 3016.7 ms. Report the achieved value, not the requested one.

Data is stored locally per machine and exported through a normal save dialog.
Both CSV shapes are long-format with participant/session/trial keys on every
row, so files `rbind` cleanly in R.

---

## When you should use something else

If you need **sub-millisecond accuracy, hardware-triggered synchronisation, or
EEG/polysomnography event markers**, use [PsychoPy](https://psychopy.org) with a
Cedrus response box or a Black Box Toolkit instead. That is the field standard
and it can do things a Chromium-based app structurally cannot.

For this paradigm the extra precision buys little: sleep latency is scored in
3-second epochs, and sleepiness-related reaction-time slowing operates on the
scale of tens to hundreds of milliseconds, comfortably above this app's residual
jitter. But if a reviewer asks for photodiode-verified timing, the honest answer
is to verify it with a photodiode — which is why the patch is built in.

---

## Reference

Bennett, L. S., Stradling, J. R., & Davies, R. J. O. (1997). A behavioural test
to assess daytime sleepiness in obstructive sleep apnoea. *Journal of Sleep
Research*, 6(2), 142–145.
