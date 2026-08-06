# BSRT — Behavioral Sleep Resistance Task

A software-based behavioral sleep resistance task, implementing the **OSLER
(Oxford Sleep Resistance Test)** paradigm as a standalone web app.

No build step, no server, no dependencies. Open `index.html` and it runs.

---

## The paradigm

OSLER is a behavioural alternative to the MWT (Maintenance of Wakefulness Test)
that requires no EEG. The participant sits in a dark, quiet room and responds to
a dim light that flashes at a fixed interval. Sleep onset is inferred
*behaviourally*, from a run of consecutive failures to respond.

Defaults reproduce the standard protocol:

| Parameter | Default | Meaning |
|---|---|---|
| Stimulus interval | 3000 ms | onset-to-onset spacing of flashes |
| Stimulus duration | 1000 ms | how long the light stays on |
| Miss criterion | 7 | consecutive missed flashes (= 21 s) scored as sleep onset |
| Maximum duration | 40 min | ceiling; trials reaching it are censored |

All four are editable on the setup screen, so the task can also be run under
modified parameters.

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

## Data

Trials are saved to **`localStorage` in that browser only**. Nothing is
transmitted anywhere. Clearing site data deletes them — **export after every
session**.

Two CSV shapes, both long-format and ready for `read.csv()` in R:

- **Epoch-level** — one row per flash: `epoch_index`, `responded`, `rt_ms`,
  `extra_responses`, plus participant/session/trial keys on every row so
  multiple files `rbind` cleanly.
- **Summary-level** — one row per trial, with all scored measures and the
  parameters the trial was run under.

"Export all" produces a single combined file across every saved trial.

### Research use

This app keeps data on the participant's device, which is privacy-friendly by
default. But once you are collecting data *from participants for research*,
consent and data-governance obligations (GDPR and your ethics approval) apply to
the exported CSVs regardless of where they were generated. This is a research
instrument, not a diagnostic device.

---

## A more accurate desktop version exists

If reaction-time precision matters to your study, use `../BSRT-desktop` instead.
It is the same paradigm packaged as a native macOS/Windows application, with
frame-locked stimulus presentation, a display calibration routine, and per-trial
timing provenance exported with the data. This browser version remains the right
choice when you need participants to run the task on their own devices without
installing anything.

## Timing accuracy

Stimuli are scheduled against **absolute targets** (`startTime + n × ISI`) rather
than by chaining relative delays, so presentation does not drift over a
40-minute trial. Reaction times are measured with `performance.now()` from the
moment the stimulus is rendered.

Browser timing is nonetheless less precise than dedicated hardware. Reaction
times carry display and input latency of a few tens of milliseconds. This is
immaterial for the primary outcome — sleep latency is scored in 3-second
epochs — but treat the RT distributions as relative measures, and do not compare
them across different devices or browsers.

The app requests a **screen wake lock** and **fullscreen** at trial start (both
are best-effort; browsers may decline). Navigating away mid-trial is counted and
reported.

---

## Running it

**Locally:** open `index.html` in a browser. That is the whole procedure.

**Hosting it** (so participants can load it on their own device): drag this
folder onto <https://app.netlify.com/drop>. It goes live at a `*.netlify.app`
URL. Sign in to claim the site — unclaimed anonymous deploys are temporary. The
same folder also works on GitHub Pages or Cloudflare Pages; all asset paths are
relative.

To update a deployed site, change the files and drag the folder onto the drop
zone on the site's Deploys page again.

---

## Reference

Bennett, L. S., Stradling, J. R., & Davies, R. J. O. (1997). A behavioural test
to assess daytime sleepiness in obstructive sleep apnoea. *Journal of Sleep
Research*, 6(2), 142–145.
