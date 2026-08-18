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
10% fastest / 10% slowest / interpercentile range for both **RT** and **RS**,
each in raw and corrected form; plus `velocity_rs`, `acceleration_rs`,
`velocity_rt`, `acceleration_rt`.

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

**Interpercentile range (IPR)** is the absolute distance between those two
decile means — `|slowest10 − fastest10|` for RT, `|fastest10 − slowest10|` for
RS. It measures within-test stability: a participant whose best and worst tenths
sit far apart is performing erratically even when the average still looks
intact. Reported per minute and for the whole test, raw and corrected
(`ipr_rt`, `ipr_rs`, `corr_ipr_rt`, `corr_ipr_rs`).

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

### Clearing the list

Two buttons sit under the picker: **Remove this participant from the list**
(the one selected, or whose ID is typed) and **Clear the whole list**. Both ask
for confirmation naming what is about to go.

**Neither deletes any data.** The roster only holds the details used to fill
the form in; recorded trials live under a separate key, stay exportable, and
are untouched. Clearing the list is also not permanent in the sense that
matters — testing someone again simply adds them back.

That first-run rebuild is a **one-time migration**, not a repair that runs
whenever the list looks empty. The difference is the whole point: refilling an
empty roster on every load would silently resurrect everyone from the stored
trials the next time the page opened, so clearing would appear to do nothing.
A flag (`…participants.v1.seeded`) records that the migration has run, and is
also set when the list is cleared on purpose, so an empty roster stays empty.

### Deleting recorded trials

Removing someone from the roster deliberately leaves their data alone. To
delete the data itself, open **Recorded trials** on the setup screen. Every
stored trial is listed with the participant, session label, trial number, when
it was run, and how it ended — enough to tell two trials of the same person
apart before removing one.

Tick the trials to delete and press **Delete selected**; the button carries the
count, and is disabled while nothing is ticked. Several trials can go at once,
and they do not have to be adjacent. A filter narrows the list to one
participant, and **Select all** then applies to what is actually shown rather
than to everything stored, so filtering to one person and selecting all cannot
sweep up somebody else's data.

Deletion asks for confirmation naming the participants affected, warns that it
cannot be undone, and suggests exporting first. **Export first.** Nothing is
recoverable afterwards: there is no undo and no copy elsewhere.

Two things deliberately do not happen. Deleting trials does not remove anyone
from the participant list — someone can be tested again without re-entering
their details, and forgetting a person is a separate decision made with the
roster buttons. And the suggested trial number is recomputed only when the
deletion touched the participant currently in the form, so removing an
unrelated trial never overwrites a number entered by hand.

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

## Device check (browser build)

The desktop build asks Electron what machine it is on. A browser cannot, so
this build measures instead — which matters most when the task is run on
someone else's laptop for screening and nobody was there to see the conditions.

Two places show it:

- **Before testing.** A **Device check** button on the setup screen measures
  the display, reads the clock resolution, and fills in the input delay as soon
  as a key is pressed. It ends with a plain verdict — *fit for testing*,
  *usable with caveats*, or *not a good device for timing-sensitive data* — and
  says which measurement drove it.
- **With every trial, automatically.** The display is calibrated during the
  three-second countdown, so it costs no extra waiting and describes the
  display in the state the task actually uses. Frames are counted throughout
  the trial, and input dispatch delay is measured from the participant's real
  responses rather than a separate key-tapping ritual.

Everything lands on the results screen under **Recording conditions**, and in
the summary CSV (`refresh_hz_measured`, `frame_interval_ms`,
`onset_quantisation_ms`, `frame_mad_ms`, `timer_resolution_ms`,
`input_dispatch_median_ms`, `dropped_rate_trial`, `ran_fullscreen`,
`device_grade`, and the rest).

### What it measures, and why each one matters

| Measure | Why it is in the file |
|---|---|
| Refresh rate, measured | Sets onset quantisation: a stimulus can only appear on a frame boundary, so 60 Hz means ±8.35 ms |
| Frame stability (MAD) | An unstable compositor makes onsets irregular even when the schedule is exact |
| Dropped frames during the trial | Reported as a rate with its denominator, since a bare count cannot be read |
| **Clock resolution** | Browsers coarsen `performance.now()` as a Spectre mitigation — 0.1 ms in Chromium, **1 ms in Firefox and Safari** unless the page is cross-origin isolated. It is a hard floor on reaction-time resolution and it is measured, not assumed |
| Input dispatch delay | How long an event took to reach the page. A Bluetooth keyboard can add tens of milliseconds. Measured and reported, **not** subtracted from reaction times |
| Ran full screen, window focus | A windowed task can be interrupted; focus lost mid-trial can look like sleepiness |
| Browser, screen, cores, memory | Provenance for comparing across devices |

### Connection latency is deliberately not a timing figure

This is the part people expect to see and it would mislead them.

**The network is not in the loop.** The whole task is client-side: once the
page has loaded, every stimulus is drawn by the local compositor and every
keypress is timed by the local clock. Nothing crosses the network between the
stimulus appearing and the response being stamped, so a participant on
satellite internet and one on fibre get identical timing. Reporting a ping
next to the reaction times would invite someone to exclude a participant for a
slow connection that had no bearing on their data.

What the connection *can* do is fail to deliver the page, deliver it slowly,
or keep loading assets while the first frames are drawn. So delivery **is**
recorded — protocol, host, load time, and the browser's own coarse connection
hint — and labelled as provenance, next to a note on screen saying exactly
this. A `file://` or `localhost` page is reported as having involved no
network at all.

### Correcting for input dispatch delay

The gap between the OS stamping a keypress and this code getting to look at it
is measured for every response, and **both readings are kept**:

| Column | What it is |
|---|---|
| `rt_event_ms` | the response measured from the OS event stamp — dispatch delay removed |
| `rt_handler_ms` | the same response from handler time — raw, delay included |
| `input_delay_ms` | their difference, for that one response |
| `rt_source_used` | which of the two was scored |

`rt_handler_ms − rt_event_ms = input_delay_ms` on every row, so the correction
can be reproduced, audited or reversed in analysis rather than taken on trust.
The results panel shows the mean of both with the difference between them and
marks the one that was scored, and the summary export carries
`rt_source_setting`, `mean_rt_event_ms`, `mean_rt_handler_ms`,
`mean_input_delay_ms`, `median_input_delay_ms`, `n_rt_compared` and
`n_stamp_rejected`.

**Device check → Time responses from** chooses which one feeds the scoring:
*OS event stamp* (default, delay removed) or *handler time* (raw). The desktop
build has the same setting and the same default, so the two builds stay
comparable.

The stamp is validated **per response**, not once per trial: a missing stamp,
or one implying a delay outside a plausible range, means this browser is not
putting `event.timeStamp` on the performance clock, so that response falls back
to handler time and is counted in `n_stamp_rejected`. Nothing is silently
substituted.

> **This changed the web build's default.** Before, reaction times were handler
> time with the delay included — the *uncorrected* variant — which meant they
> were not directly comparable with the desktop's. They now default to the same
> corrected reading. On a wired keyboard the shift is well under a millisecond;
> on a Bluetooth one it can be tens. Set **handler time** to reproduce the old
> behaviour, and note that `rt_handler_ms` reproduces it in any case.

### One thing worth knowing about the frame monitor

The browser build schedules stimuli with `setTimeout`, and the frame monitor
only observes. Measured over repeated runs it does not degrade onset accuracy
— it slightly *improves* it, because a live `requestAnimationFrame` loop stops
Chromium coarsening timers on a page it thinks is idle. With the monitor
running, the worst onset error over a trial was 0.5–0.8 ms; with it disabled,
the same test drifted to 3.9 and 11.4 ms.

---

## Preliminary norms

**Off by default.** The results screen shows the plain, uncoloured tables it
always did — plus the new interpercentile range rows. A single switch,
**Show preliminary norms**, reveals the comparison, and the choice is
remembered so it does not have to be re-ticked every trial. The comparison is
always computed and always exported; only the display is optional.

When switched on, results are compared against a **preliminary** reference of
435 test sessions from 48 people, pooled from **two studies**, and binned by the
**hour the session started** and by **cumulative test length**. A 3-minute test
is compared against the first 3 minutes of the reference sessions, never against
their whole run.

| Study | Sessions | People | Protocol | Hours covered |
|---|---|---|---|---|
| Control group | 291 | 12 | 8 minutes | all 24 |
| Overnight study | 144 | 36 | 40 minutes | 23:00, 01:00, 03:00, 05:00 |

So the four overnight hours pool both studies for lengths 1–8 minutes and carry
the overnight study alone from 9 to 40; **every other hour stops at 8 minutes**.
The app picks the deepest window the hour actually supports, so a 20-minute test
at 14:00 is compared over its first 8 minutes, and the same test at 03:00 over
all 20.

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

### The two samples are not interchangeable

This matters more than any other caveat here, and the app states it on screen
whenever a pooled or overnight-only cell is used.

Over the same four hours, the overnight sample is **slightly faster** than the
control group (RT average 316 ms vs 333 ms at 8 minutes) but misses **about nine
times as many trials** (2.85% vs 0.32%), and **52 of its 144 sessions (36%)**
reached the seven-consecutive-miss sleep-onset criterion and stopped early — a
thing that essentially never happened in the control group. Fast when awake,
frequently absent altogether, is the signature of a sleepy population rather
than a merely tired one.

That shape is worth confirming before these numbers are quoted as norms. If the
overnight sample is a clinical group, hours 23/01/03/05 are not a healthy
reference and should be split back out. Nothing is lost either way: every cell
records its own composition as `prov[hour][length-1]` = `[control sessions,
overnight sessions, sessions that had already ended at sleep onset]`, those
counts are exported with every comparison row, and the pooling can be reversed
by regenerating from the raw workbooks.

### Censoring at the longer lengths

A session contributes to every length it completed and to none after it. Because
sessions that reach sleep onset stop, they drop out of the longer windows:

| Window | Sessions contributing | Already ended at sleep onset |
|---|---|---|
| 8 min | 429 | 6 |
| 10 min | 133 | 11 |
| 20 min | 117 | 27 |
| 30 min | 99 | 45 |
| 40 min | 92 | 52 |

By 40 minutes the reference describes only the people who stayed awake. That is
a **survivor bias in the direction of alertness**, and it runs the wrong way for
screening: it makes a sleepy participant look worse than the true population
would. The app names the number of dropped-out sessions whenever it uses such a
cell, and the count is in the export.

### When no comparison is offered

The reference sessions are BSRT runs at a 3000 ms interval, which is 20 stimuli
a minute. The panel refuses, and says which condition failed, when the trial is
a **PVT** (about ten stimuli a minute), when the **interval is not 3000 ms**,
when the **start hour** cannot be read, or when the test did not complete a
**full minute**. Trial, miss and lapse counts would otherwise be compared
against a reference measuring something else.

### Limits you should read before quoting a z-score

- **They are preliminary.** Most hour bins hold only **9–15 sessions**, and the
  same participant contributes to many bins. The observations are not
  independent, and every SD is itself estimated from a small number of values.
- **Two samples, two expected durations.** The control participants knew they
  were settling in for 8 minutes and the overnight participants for 40. Someone
  expecting a short test may pace themselves differently, so **short comparisons
  can be biased for that reason alone**. The app flags this whenever the window
  is shorter than a contributing study's protocol.
- **Some reference cells have almost no spread.** Every control session scored
  20 hits out of 20 in the first minute, so a difference of one or two trials
  can produce a double-digit SD figure. Where the reference SD is exactly zero
  no z-score is computed at all; the value is marked as outside the reference
  range instead of dividing by zero.
- **The comparison uses the reference conventions, not the app's own scoring.**
  Every response counts however slow, and a decile is sized with `ceil(10%)`
  rather than `round(10%)`. **The numbers in the norms panel will not always
  match the tables above it**, which is deliberate — the export keeps both.

### Verification

The tables are computed from the **raw trial values** of both studies, not from
anyone's published means, by a Python port of the same `normativeSummary()` the
app applies to a participant. Three checks back that:

- Recomputing the previous 24 × 8 published tables from the control group's raw
  trials reproduced **all 8064 cells to within 5 × 10⁻⁵**, which is the rounding
  of the published file. The earlier norms were correct.
- Each workbook's own derived columns (`Mean_*`, `SD_*`, `Md_*`, `Miss_*`,
  `FS_*`) were reproduced from its raw trials — **291/291** control sessions and
  **all six test lengths** of the overnight study — which fixes the conventions
  rather than assuming them. Note those columns keep timeouts in the RT average;
  the app's norms exclude them, as the app's own scoring does.
- Every populated cell is checked to be 42 values wide with `n` equal to the sum
  of its provenance counts: 320 cells, no exceptions.

### Export

A fourth CSV, **norms comparison**, is long-format with one row per variable:
`norm_value`, `norm_ref_mean`, `norm_ref_sd`, `z_worse`, `band`, plus the hour,
window, reference n, the `norm_below_protocol` / `norm_truncated` flags, and the
cell's composition as `norm_ref_n_control`, `norm_ref_n_overnight`,
`norm_ref_n_ended_early`, `norm_mixed_studies` and `norm_hour_max_min`. It
`rbind`s with the other exports.

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

## There is also a PsychoPy version

`../BSRT-psychopy` runs the same task as a PsychoPy experiment, for labs that
already have PsychoPy set up. It takes stimulus onsets from the vsync timestamp
of the frame that carried them and reaction times from the keyboard's own
hardware stamp, which is more than either of the other builds can do.

Same paradigm, same scoring, same CSV columns — the Python scoring is a checked
port of `scoring.js` and a test drives both with identical inputs and asserts
they agree to 1e-9, so the three builds' data can be pooled.

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

**Getting a copy without git:** download `downloads/BSRT-webapp.zip` from the
repository — click the file, then the Download button — and unzip it. That is
this folder, bundled for people who should not have to use git.

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
