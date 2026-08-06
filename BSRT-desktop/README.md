# BSRT (desktop) — Behavioral Sleep Resistance Task

An OSLER-paradigm sleep resistance task packaged as a native application for
**macOS and Windows**, with frame-locked stimulus presentation and a calibration
routine that measures and reports the timing characteristics of the machine it
is running on.

This is the accuracy-focused sibling of the browser version in `../BSRT`. The
paradigm, scoring, and CSV shape are the same. What differs is how time is
measured.

---

## Building the executables

You need [Node.js](https://nodejs.org) 18 or newer. Everything else is fetched
by `npm install`.

```bash
cd BSRT-desktop
npm install
npm start          # run it without packaging, to check it works
```

Then build for your platform:

```bash
npm run dist:mac   # -> dist/BSRT-1.0.0-arm64.dmg, -x64.dmg, and .zip
npm run dist:win   # -> dist/BSRT Setup 1.0.0.exe and a portable .exe
```

**You must build macOS apps on a Mac and Windows apps on Windows.** Code signing
and platform packaging are not portable. If you need both from one machine, run
the builds in CI (GitHub Actions with a `macos-latest` and a `windows-latest`
job) rather than trying to cross-compile.

### The unsigned-app warnings

These builds are not code-signed, because signing requires a paid Apple
Developer account and a Windows code-signing certificate. Participants will
therefore see a scary dialog on first launch:

- **macOS:** right-click the app → *Open* → *Open*. Or once, in Terminal:
  `xattr -cr /Applications/BSRT.app`
- **Windows:** SmartScreen shows "Windows protected your PC" → *More info* →
  *Run anyway*.

If you are distributing to participants rather than running the task yourself,
budget for signing certificates — talking people through a security warning over
the phone is a poor start to a testing session.

---

## Languages

The interface runs in **English, French, Dutch and German**, switchable from the
selector at the top. The language a participant saw is exported as `language`
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

When the sleep-onset criterion is reached, the app sounds a three-tone alarm.
It is on by default and can be switched off; the setting is exported as
`alarm_enabled`. The tone is synthesised rather than played from a file, so the
app stays dependency-free and works offline.

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
`epoch_index`, `minute`, `onset_ms`, `responded`, `rt_ms`, `rs_per_sec`,
`outcome`, `lapse`, `late_response`, `anticipation`, `extra_responses`.

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

Name, address and birth date are recorded and written into every export. They
are directly identifying. For research use, prefer a pseudonymous participant ID
and hold the identity key separately under your data-protection plan — the app
warns about this on the setup screen. Date and time of testing are captured
automatically.

### Scoring code

All derived metrics come from `scoring.js`, which is byte-identical between the
browser and desktop builds so the two cannot diverge.
`npm run check:scoring` in `BSRT-desktop` asserts the copies match.

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

Each interval is rounded independently to a whole number of frames, and the
boundaries accumulate, so a variable schedule is delivered as precisely as a
fixed one — measured deviation was under 0.1 ms per interval in testing.

Because scheduling is by frame count, a **dropped frame shifts everything after
it** by one frame interval in wall-clock terms while keeping the stimulus
periodic on screen. That is the deliberate trade: the participant sees a regular
rhythm, and the drift is bounded by the dropped-frame count, which is measured
and exported (`dropped_frames_trial`). On a healthy display it is zero. If it is
not, the display check will already have graded the machine `poor`.

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
